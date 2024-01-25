%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulates Selection Gene Drive System %%%
%%%         Scott Leighow 09/22/21        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all
% rng(3)

% 09/16/22 Update: added birth and death inputs for sensitivity analysis

function output = SRG_20240110(pop_treat,mut,inf_eff,swtch_dly,dim_res,a_B,b,d,npats)

%% Parameters

%%% Rates [/day] %%%
% (Bozic eLife 2014)

% b = 0.14;   % division rate
% d = 0.13;   % natural death rate
a_T = 0.04; % targeted therapy drug kill rate
a_X = 0.0; % direct toxic metabolite drug kill rate (not considering bystander effect; only affects cells expressing suicide gene)
% a_B = 0.03; % toxic metabolite bystander effect drug kill rate (effect on all cells sensitive to metabolite, not just those expressing suicide gene)
% *Note that a_B will ultimately be proportional to the fraction of tumor
% cells expressing suicide gene

% mut = 1e-8;     % Mutation rate [/bp/division]

%%% Populations %%%

% P = [S R M D G Q L C]

% S = sensitive
% R = resistant to targeted therapy
% M = resistant to toxic metabolite
% D = resistant to both
% G = infected
% Q = infected and resistant to targeted therapy
% L = infected and resistant to toxic metabolite
% C = infected and resistant to both

%%% Initial conditions %%%

S_0 = 100;

pop_init = [S_0 0 0 0 0 0 0 0];
npop = length(pop_init);

%%% Simulation conditions %%%

% inf_eff = 0.1;      % infection efficiency - what percentage of tumor cells take up gene drive?
% pop_treat = 1e9; % size of population upon treatment (also determines tumor size upon switch activation and relapse)
% swtch_dly = 360;    % delay between start of switch 2 and end of switch 1

tau_init = NaN;     % length of model time intervals [days] (if NaN model determines tau adaptively)
epsilon = 0.03;     % error tolerance []
% npats = 1;

%% Simulations

%%% Mutation matrix 
%%% Each element represents number of sites for mutation (from = row; to = cols)

%          S R M D G Q L C
mut_mat = [0 1 1 0 0 0 0 0   % S
           0 0 0 1 0 0 0 0   % R
           0 0 0 1 0 0 0 0   % M
           0 0 0 0 0 0 0 0   % D
           1 0 0 0 0 1 1 0   % G
           0 1 0 0 0 0 0 1   % Q
           0 0 1 0 0 0 0 1   % L
           0 0 0 1 0 0 0 0]; % C
       
%%% Sensitivity vectors
%%% Which populations are sensitive to which therapies?

% To what degree does dimerizer confer resistance to gene drive cells?
% dim_res is a scaling factor for the targeted therapy drug kill rate in
% infected cells
% 0=completely resistant; 1=completely sensitive
% neutral ngr at (b-d-a_T*dim_res)=0 ==> dim_res=(b-d)/a_t=0.25 for above
% parameters
% dim_res = 0;  
       
% Targeted therapy
%            S R M D G Q L C
Tsen_vec0 = [0 0 0 0 0 0 0 0];  % sensitivity to targeted therapy (before switch 1 - so no effect)
Tsen_vec1 = [1 0 1 0 dim_res 0 dim_res 0];  % sensitivity to targeted therapy (switch 1); 1 = sensitive, 0 = resistant
Tsen_vec2 = [1 0 1 0 1 0 1 0];  % sensitivity to targeted therapy (switch 2)

% Toxic metabolite direct killing (without bystander effect, i.e. only 
% cells that express suicide gene)
Xsen_vec0 = [0 0 0 0 0 0 0 0];  % sensitivity to direct toxic metabolite (before switch 2 - no effect)
Xsen_vec2 = [0 0 0 0 1 1 0 0];  % sensitivity to direct toxic metabolite (switch 2)

% Toxic metabolite bystander killing (i.e. all cells regardless of suicide
% gene expression)
Bsen_vec0 = [0 0 0 0 0 0 0 0];  % sensitivity to direct toxic metabolite, not bystander effect (before switch 2 - no effect)
Bsen_vec2 = [1 1 0 0 1 1 0 0];  % sensitivity to toxic metabolite (switch 2)

% Population groups (needed for some calculations, e.g. strength of
% bystander effect)
sg_pops = [0 0 0 0 1 1 1 1];    % expresses suicide gene
res_pops = logical([0 1 1 1 0 1 1 1]);   % which populations are generally resistant? monitor for relapse

   
%%% Initialize storage variables
output = NaN(npats,2); % columns = [eradicated? PFS/time2erad]

% timeSet = cell(npats,1);
% cellSet = cell(npats,1);

parfor r = 1:npats
        
    % Initialize
    n = 1;
    t = 0;
    tau = tau_init;
    singleSSAs = NaN;
    pop = pop_init;
    
    % Define switch states and associated drug kill rates
    switch1 = false;
    switch2 = false;
    t_switch1end = NaN;
    Tsen = Tsen_vec0;
    Xsen = Xsen_vec0;
    Bsen = Bsen_vec0;
    
    %%% Determine Jacobian and numChange (reflects propensity vector below)
    %%% Note: can be determined outside of simulation, since they are 
    %%% (mostly) state-independent

    % Determine Jacobian for mutational events  
    Jac_mut = zeros(npop*npop,npop);
    for i = 1:npop
        rows_i = (8*i-7):8*i;
        Jac_mut(rows_i,i) = mut*b*mut_mat(i,:)';
    end

    % Compile Jacobian
    Jac = [b*eye(npop);         % division
           Jac_mut;             % mutation
           d*eye(npop);         % natural death
           diag(a_T*Tsen);      % targeted therapy killing
           diag(a_X*Xsen);      % toxic metabolite direct killing
           diag(0.5*a_B*Bsen)]; % toxic metabolite bystander killing - note: rather than calculate derivative of suicide gene fraction term, recognize that fraction is between 0 and 1, so assume 0.5 for sake of computational efficiency

    numChange = [eye(npop);                     % faithful division
                 repmat(eye(npop),npop,1);      % mutation
                 repmat(-eye(npop),npop,1)];    % all death events

    % Run simulation until resistant population reaches detection (relapse)
    while sum(pop(end,:))>0 && sum(pop(end,res_pops))<pop_treat
        
        % Define current
        P_curr = pop(end,:);        
        t_curr = sum(t);

        %%% Treatment conditions
        
        % Switch 1
        if sum(P_curr) >= pop_treat && ~switch1
            
            switch1 = true;
            
            % Infect cells
            pop(end,5:8) = round(pop(end,1:4)*inf_eff);
            P_curr = pop(end,:);
            
            % Begin treating with targeted therapy and dimerizer
            Tsen = Tsen_vec1;
            
            % Recalculate Jacobian
            Jac = [b*eye(npop);         % division
                   Jac_mut;             % mutation
                   d*eye(npop);         % natural death
                   diag(a_T*Tsen);      % targeted therapy killing
                   diag(a_X*Xsen);      % toxic metabolite direct killing
                   diag(0.5*a_B*Bsen)]; % toxic metabolite bystander killing - note: rather than calculate derivative of suicide gene fraction term, recognize that fraction is between 0 and 1, so assume 0.5 for sake of computational efficiency

               
            t_treat = t_curr;
            
        end
        
        % Switch1/2
        if P_curr(5) >= pop_treat && switch1 && ~switch2
            
            switch2 = true;
            
            % Begin treating with prodrug
            switch2 = true;
            Xsen = Xsen_vec2;
            Bsen = Bsen_vec2;
            
            % Recalculate Jacobian
            Jac = [b*eye(npop);         % division
                   Jac_mut;             % mutation
                   d*eye(npop);         % natural death
                   diag(a_T*Tsen);      % targeted therapy killing
                   diag(a_X*Xsen);      % toxic metabolite direct killing
                   diag(0.5*a_B*Bsen)]; % toxic metabolite bystander killing - note: rather than calculate derivative of suicide gene fraction term, recognize that fraction is between 0 and 1, so assume 0.5 for sake of computational efficiency

            t_switch1end = t_curr+swtch_dly;
            
        end    
        
        % Switch 2
        if t_curr>=t_switch1end
            
            % Remove dimerizer
            Tsen = Tsen_vec2;
            
            % Recalculate Jacobian
            Jac = [b*eye(npop);         % division
                   Jac_mut;             % mutation
                   d*eye(npop);         % natural death
                   diag(a_T*Tsen);      % targeted therapy killing
                   diag(a_X*Xsen);      % toxic metabolite direct killing
                   diag(0.5*a_B*Bsen)]; % toxic metabolite bystander killing - note: rather than calculate derivative of suicide gene fraction term, recognize that fraction is between 0 and 1, so assume 0.5 for sake of computational efficiency

        end
        
        %%% Calculate Propensitiy vector
        
        % Calculate probabilities of mutation events
        mut_vec = reshape(mut_mat',[],1); % rearrange mut_vec to match dimensions of evts variable (n x 1)
        mut_evts = b*mut*mut_vec.*repelem(P_curr',npop); % gives probability of each mutation event in vector format; begins with all mutations from S, then R, etc.
        
        % Calculate probabilities of targeted therapy kill events
        T_evts = a_T*Tsen'.*P_curr';
        
        % Calculate probabilities of toxic metabolite direct kill events
        % (no bystander effect)
        X_evts = a_X*Xsen'.*P_curr';
        
        % Caclulate probabilities of toxic metabolite bystander kill events
        % (proportional to fraction of tumor expressing suicide gene)
        bystndr = (sg_pops*P_curr')/sum(P_curr); % what fraction of cells express suicide gene?
        B_evts = a_B*bystndr*Bsen'.*P_curr';
        
        % Define propensity vector and change matrix 
        evts = [b*P_curr';  % faithful division
                mut_evts;   % mutation
                d*P_curr';  % natural death  
                T_evts;     % targeted therapy killing
                X_evts;     % toxic metabolite direct killing
                B_evts];    % toxic metabolite bystander killing

        % Calculate total rate
        theta = sum(evts);
                        
        % Tau-leaping model
        if isnan(singleSSAs)

            % Determine tau if not defined
            if isnan(tau)

                % Calculate tau
                fjj = Jac*numChange';
                mu = fjj'*evts;
                std = (fjj.^2)'*evts;

                tau = min(min(epsilon*theta./(abs(mu)), epsilon^2*theta^2./(std)));
            end

            % Rerun with singleSSAs if tau is too small
            if tau < 1/theta
                singleSSAs = 10;
                tau = tau_init;
                continue

            % Otherwise calculate new state
            else
                % For each event, find number of times it's occured
                % in time interval tau
                P_next = P_curr;
                for i = 1:length(evts)
                    numEvt = poissrnd(evts(i)*tau);
                    P_next = P_next + numChange(i,:).*numEvt;
                end
                P_next = round(P_next);

                % Rerun with smaller tau if any populations are
                % negative
                if sum(P_next < 0) > 0
                    tau = tau/2;
                    continue

                % Otherwise save results in array
                else                      
                    timeElapsed = tau;
                    t(n+1) = timeElapsed;
                    pop(n+1,:) = P_next;

                    % Reset tau
                    tau = tau_init;
                end
            end

        % Gillespie without tau leaping
        else
            r1 = rand;
            timeElapsed = -log(r1)/theta;
            t_curr = sum(t(1:n));

            t(n+1) = timeElapsed;

            idx = randsample(length(evts), 1, true, evts);
            P_next = round(P_curr + numChange(idx,:));
            pop(n+1,:) = P_next;

            % Reduce singleSSAs
            singleSSAs = singleSSAs - 1;
            if singleSSAs < 1
                singleSSAs = NaN;
            end
        end
                        
        n = n + 1;
    end
    
    % In the case of eradication, repeat if occured pretreatment
    if sum(pop(end,:)) == 0 && ~switch1
        continue
    end
    
    % Save time vector as cumulative sum of time intervals
    time = cumsum(t)';
    
    % Save results
    
    erad = sum(pop(end,:))==0;
    t2end = t_curr - t_treat;
    
    output(r,:) = [erad t2end];
    
    disp(r);
            
end

end



%% Plot Results
% 
% colors = [0 0 1; rand(npop-1,3)];
% 
% figure
% 
% for j = 1:npop
%     semilogy(time/365,pop(:,j),'Color',colors(j,:),'LineWidth',2);
%     hold on
% end
% 
% title('Pop Dynamics')
% xlabel('Time [years]')
% ylabel('Population Size [cells]')
% legend({'S' 'R' 'M' 'D' 'G' 'Q' 'L' 'C'})
% hold off
% % 
% % subplot(1,2,2);
% for i = 1:npats
%     for j = 1:M
%         LSCs = cellSet{i}(:,1:M);
%         WBCs = cellSet{i}(:,3*M+1:4*M);
%         semilogy(timeSet{i}/365,LSCs(:,j),'o','Color',colors(j,:))
%         hold on
%     end
% end
% title('White Blood Cell Dynamics')
% xlabel('Time [years]')
% ylabel('Population Size [cells]')
% hold off
% 
% figure
% for j = 1:Rsites
%     plot(biasR(j),IC50R(j),'.','Color',colors(j+1,:),'markers',25)
%     hold on
% end
% 
% title('Parameter Space')
% xlabel('Probability')
% ylabel('IC_{50} [nM]')
% hold off
% 
% max(cellSet{i}(:,1))
