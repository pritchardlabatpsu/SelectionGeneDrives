%%% Spatial Agent-Based Model of Gene Drive System
%%% Scott Leighow - 07/09/22

% 09/02/22 Update: Some (actually quite a few) iterations in the last
% simulation were reaching the simulation tmax before eradication, even
% though it looked like a lot of them were gonna get there. So updating so
% that tmax is now 750 (instead of 365) and the simulation ends when either
% a) all cells are gone
% b) 30 days have passed since gene drive cells are eliminated (to see if
% resistance has started to grow out; gddepleted and gddepleted30
% variables)
% c) when tmax is reached
% d) when nmax is reached (now 1.3*n0 rather than 2*n0)
% Essentially, we want to extend the tmax so that there's still a limit on
% the simulation in case things go wrong, but given how long a simulation
% can run now, we want to make sure it only runs for as long as it needs to
% (i.e. until it's obvious whether or not eradication was achieved).

% 09/03/22 Update: Increased gddepleted time from 30 to 150. Increased tmax
% from 700 to 1000. Make swtch12_delay an input variable.

% 09/03/22c Update: Fixed bug in code. Previous version had error when
% determining adj_mat (roughly lines 242 and 277), where it would grab the
% kth index of the pos_curr, when it should be grabbing gdv_idx(k)th index.

% clear
% close all
% rng(30)

function [output] = GeneDriveABMfunc_090322c(size0,tight,n_foci,kill_rad,swtch12_delay)

a_B = 0.03; % toxic metabolite bystander effect
% size0 = 3000; % initial population size
res_frac = 0.005; % resistant population frequency
gd_frac = 0.05; % gene drive population frequency
% kill_rad = 5; % radius of killing for gene drive cells 
% tight = 0.3; % measure of how closely spaced gene drive cells are; large numbers give a very tightly-spaced focus of infecion; 0 gives random spatial distribution
% n_foci = 3; % number of foci of gene drive infections
n_swtch2 = 0.8*size0; % threshold number of gene drive cells to trigger switch 2
bdgmax = Inf; % maximum budging distance; defines constraint of tumor on cell division
clrnce = 0.3; % clearance rate of activated metabolite [/day]; determines how long bystander effect remains after death of gene drive cell
% swtch12_delay = 180; % delay between turning on switch 2 and turning off switch 1 (protects gene drive cells from TKI killing for a period of time)

%% Parameters

% Rates [/day]
% Bozic et al eLife 2013
b = 0.14;   % division rate
d = 0.13;   % natural death rate
a_T = 0.04; % targeted therapy drug kill rate
a_X = 0.0; % direct toxic metabolite drug kill rate (not considering bystander effect; only affects cells expressing suicide gene)
% a_B = 0.03; % toxic metabolite bystander effect drug kill rate (effect on all cells sensitive to metabolite, not just those expressing suicide gene)
% *Note that a_B will ultimately be proportional to the fraction of tumor
% cells expressing suicide gene
                                
% Population parameters
nsen0 = round(size0*(1-res_frac-gd_frac)); % number of sensitive cells
nres0 = round(size0*res_frac);  % number of resistant cells
ngdv0 = round(size0*gd_frac);  % number of gene drive cells

if n_foci > ngdv0
    warning('Number of infection foci exceeds number of infected cells')
end

n0 = sum([nsen0 nres0 ngdv0]);

% Define neighborhood matrix
nghbrs1 = combvec(-1:1,-1:1,-1:1)';
% [~,row0s] = ismember([0 0 0],nghbrs1,'rows');
% nghbrs1(row0s,:) = [];

% Simulation parameters
nmax = round(1.3*n0); % maximum number of cells

% output = NaN(nsim,4);

% Pretty colors
clr_sen = [42 186 252]/255;
clr_res = [238 29 36]/255;
clr_gdv = [26 211 26]/255;
% clr_ads = [143 100 255]/255;
% clr_adr = [143 0 255]/255;
clr_ads = clr_sen;
clr_adr = clr_res;

%% Initialize Tumor
    
% Randomly seed cells
len = ceil(n0^(1/3)); % determine length of initial cube
len = len + ceil(log10(size0)) + 1; % Add additional space based on size0 to make geometrically realistic
if rem(len,2)==0       % if length is even, increase by one
    len = len+1;
end
lenvec = -(len-1)/2:(len-1)/2;
cube0 = combvec(lenvec,lenvec,lenvec)';    
unocc_idx = (1:size(cube0,1))';

dst2foci = sqrt(sum(cube0.^2,2)); % calculate Euclidean distance of potential points
pos0 = [0 0 0];
while size(pos0,1) < n0 % seed cells (preferentially close to center) until all cells are assigned position
    unocc_idx = (1:size(cube0,1))';
    unocc_idx(ismember(cube0,pos0,'rows')) = [];
    cells_idxi = randsample(unocc_idx,n0-size(pos0,1),true,exp(-5*dst2foci(unocc_idx)));
    pos_i = cube0(cells_idxi,:);
    pos0 = [pos0; pos_i];
    pos0 = unique(pos0,'rows');
end

% Distribute gene drive cells
% Gene drive cells are dispersed by first choosing a randomly assigned
% center
mean_pos_gdv = pos0(randsample(1:n0,n_foci),:);

% Cells are "infected" randomly, with weight given to cells closest to
% center of infection (exponential distribution) - similar to assignment of
% cell positions above
dst2foci = NaN(n0,n_foci);
for n = 1:n_foci
    dst2foci(:,n) = sqrt(sum((pos0-mean_pos_gdv(n,:)).^2,2)); % calculate Euclidean distance of cells from infection foci
end
dst2focus = min(dst2foci,[],2); % select distance from nearest foci

gdv0 = mean_pos_gdv;

while size(gdv0,1) < ngdv0 % seed cells until all gene drive cells are assigned a position
    uninfctd_idx = 1:n0;
    uninfctd_idx(ismember(pos0,gdv0,'rows')) = [];
    cells_idxi = randsample(uninfctd_idx,ngdv0-size(gdv0,1),true,exp(-tight*dst2focus(uninfctd_idx)));
    gdv_i = pos0(cells_idxi,:);
    gdv0 = [gdv0; gdv_i];
    gdv0 = unique(gdv0,'rows');
end

% Establish phenotypes:
% 1 = sensitive (sen)
% 2 = resistant (res)
% 3 = gene drive (gdv)
% 4 = sensitive, adjacent to gene drive (ads)
% 5 = resistant, adjacent to gene drive (adr)

phn0 = ones(n0,1);
phn0(ismember(pos0,gdv0,'rows')) = 3;

% Distribute resistant cells (randomly dispersed)

uninfctd_idx = find(phn0==1);
res_idx = randsample(uninfctd_idx,nres0,false);
phn0(res_idx) = 2;

% Define gene drive kill radius (determines adjacent population)
kill_rad_mat = combvec(-kill_rad:kill_rad,-kill_rad:kill_rad,-kill_rad:kill_rad)'; % create matrix of neighbors up to diff_rad units away
kill_rad_mat(sqrt(sum(kill_rad_mat.^2,2))>kill_rad,:) = []; % remove positions more than diff_rad units away

% % Visualize initial spheroid
% clrs = NaN(n0,3);
% clrs(phn0==1,:) = repmat(clr_sen,sum(phn0==1),1);
% clrs(phn0==2,:) = repmat(clr_res,sum(phn0==2),1);
% clrs(phn0==3,:) = repmat(clr_gdv,sum(phn0==3),1);
% clrs(phn0==4,:) = repmat(clr_ads,sum(phn0==4),1);
% clrs(phn0==5,:) = repmat(clr_adr,sum(phn0==5),1);
% scatter3(pos0(:,1),pos0(:,2),pos0(:,3),5e2,clrs,'.')
% hold on
% scatter3(mean_pos_gdv(1),mean_pos_gdv(2),mean_pos_gdv(3),5e2,'k*')
% hold off

%% Initialize Simulation

% Initialize state
pos_curr = pos0;
phn_curr = phn0;
state_curr = [pos0 phn0];

n_curr = n0;
pop0 = [sum(phn0==1)
        sum(phn0==2)
        sum(phn0==3)
        sum(phn0==4)
        sum(phn0==5)]';
n_rsdl = 0;

% Initialize treatment conditions
alpha_T = [a_T 0 0 a_T 0]; % under switch 1, sen and ads are sensitive to TKI
alpha_X = zeros(1,5);      % under switch 1, no fitness cost of switch 2 production
alpha_B = zeros(1,5);      % under switch 1, no bystander effect

swtch1 = true;
swtch2 = false;
t_swtch2 = Inf;

bystndr_intrvl = round(size0); % bystndr_mat is updated every bystndr_intrvl iterations - general rule of thumb: every size0 iterations works well

% Initialize simulation variables
j = 1;
t = 0;
t_curr = sum(t);
tmax = 1000;

gddepleted = false;     % used to check if gene drive population has been eliminated
gddepleted150 = false;   % used to check if 10 days have passed since gene drive population has been depleted (condition in while loop)

% Initialize storage matrices
% state = NaN(nmax,4);
% state(1:size(state_curr,1),:) = state_curr;
s = 1;

output = [j t_curr pop0];
if size0 >= 1e3
    str_intrvl = round(size0/20);   % store every x iterations in output matrix
else
    str_intrvl = 5;
end

%% Simulation

%     while size(state_curr,1) > 1 && size(state_curr,1) < nmax && sum(t) < 30
% while sum(phn_curr==2)>0 && sum(phn_curr==2)<size0 && size(state_curr,1) < nmax
while n_curr>0 && n_curr<nmax && t_curr<tmax && ~gddepleted150

    % Engage switch 2 next iteration when condition met
    if sum(phn_curr==3) > n_swtch2 && ~swtch2
        swtch2 = true;
        t_swtch2 = t_curr;
        
        % Update kill rates
        alpha_X = [0 0 a_X 0 0];     % under switch 2, fitness cost of switch 2 production
        alpha_B = [0 0 a_B a_B a_B]; % under switch 2, bystander effect (only affects gdv, ads, and adr cells)
        % alpha_T updated after switch 1 turned off

        % Determine bystndr_mat
        % The variable bystndr_mat stores the positions of all sites within
        % the kill radius of a gene drive cell (first 3 cols) and the
        % number of gene drive cells within a kill radius (4th column).
        % Every time a gene drive cell dies, 1 is subtracted from the 4th
        % column of all sites in bystndr_mat that it was associated with.
        % Sites that get down to 0 are cleared at a rate of clrnce. Given
        % the budging of cells, bystndr_mat may not be exact, but this
        % saves from having to calculate the number of adjacent sites fresh
        % every iteration. A new, updated bystndr_mat is calculated every
        % bystndr_intrvl iterations.
        adj_mat = [];
        gdv_idx = find(phn_curr==3); % index of all gene drive cells
        for k = 1:length(gdv_idx)
            pos_gdvk = pos_curr(gdv_idx(k),:); % pick kth gene drive cell
            adj_k = pos_gdvk+kill_rad_mat; % identify all positions within kill radius of kth gene drive cell
            adj_mat = [adj_mat;adj_k]; % add positions to bystndr_mat
        end

        % Collapse bystndr_mat to a matrix with unique rows and the number
        % of gene drive cells associated with each row (number of times
        % that unique row appeared in the uncollapsed matrix).
        [ii,jj,kk] = unique(adj_mat,'rows','stable');
        freq = histc(kk,1:numel(jj)); % determine frequency of each row
        bystndr_mat = [ii freq];      % save in bystndr_mat

        % Update phenotype for cells within kill_rad of gene drive cells
        phn_curr(ismember(pos_curr,bystndr_mat(:,1:3),'rows')&phn_curr==1) = 4; % assign all sen in bystndr_mat to ads
        phn_curr(ismember(pos_curr,bystndr_mat(:,1:3),'rows')&phn_curr==2) = 5; % assign all res in bystndr_mat to adr
        
    end

    % Determine when to turn off switch 1 (possibility for delay after
    % engaging switch 2)
    if t_curr > t_swtch2+swtch12_delay && swtch1
        swtch1 = false;
        % Update kill rates
        alpha_T = [a_T 0 a_T a_T 0]; % after turning off switch 1, sen, gdv ads are sensitive to TKI
    end
    
    % Every bystndr_intrvl iterations, calculate fresh bystndr_mat
    if swtch2 && rem(j,bystndr_intrvl) == 1
        
        rsdl_bystndr = bystndr_mat(bystndr_mat(:,4)==0,:); % note residual rows of bystndr_mat (those yet to be cleared)

        % Determine adj_mat
        adj_mat = [];
        gdv_idx = find(phn_curr==3); % index of all gene drive cells
        for k = 1:length(gdv_idx)
            pos_gdvk = pos_curr(gdv_idx(k),:); % pick kth gene drive cell
            adj_k = pos_gdvk+kill_rad_mat; % identify all positions within kill radius of kth gene drive cell
            adj_mat = [adj_mat;adj_k]; % add positions to bystndr_mat
        end
        
        % Collapse bystndr_mat to a matrix with unique rows and the number
        % of gene drive cells associated with each row (number of times
        % that unique row appeared in the uncollapsed matrix).
        [ii,jj,kk] = unique(adj_mat,'rows','stable');
        freq = histc(kk,1:numel(jj)); % determine frequency of each row
        bystndr_mat = [ii freq;       % save in bystndr_mat
                       rsdl_bystndr]; % include rows from previous iteration of bystndr_mat yet to be cleared (residual rows)

        % Update phenotypes
        phn_curr(phn_curr==4) = 1; % reassign all ads to sen
        phn_curr(phn_curr==5) = 2; % reassign all adr to res
    
        if size(bystndr_mat,1)>0
            phn_curr(ismember(pos_curr,bystndr_mat(:,1:3),'rows')&phn_curr==1) = 4; % assign all sen in bystndr_mat to ads
            phn_curr(ismember(pos_curr,bystndr_mat(:,1:3),'rows')&phn_curr==2) = 5; % assign all res in bystndr_mat to adr
        end

    end
    
    % Current populations
    nS = sum(phn_curr==1);
    nR = sum(phn_curr==2);
    nG = sum(phn_curr==3);
    nAS = sum(phn_curr==4);
    nAR = sum(phn_curr==5);
    pop_curr = [nS nR nG nAS nAR];

    if swtch2 && size(bystndr_mat,1)>0
        n_rsdl = sum(bystndr_mat(:,4)==0); % number of positions with residual bystander effects (but without a living gene drive cell)
    else
        n_rsdl = 0;
    end
        
    % Propensity vector
    evts = [(b.*pop_curr)';         % cell division event
            (d.*pop_curr)';         % natural death event
            (alpha_T.*pop_curr)';   % TKI kill event
            (alpha_X.*pop_curr)';   % direct toxic metabolite kill event (not bystander effect; fitness cost of switch 2 production)
            (alpha_B.*pop_curr)';   % bystander effect kill event
            n_rsdl*clrnce];         % bystander effect clearance event (remove residual bystander effect) 

    % Draw event and time to event
    theta = sum(evts);
    t(j) = -log(rand)/theta;

    idx_evt = randsample(length(evts),1,true,evts);

    % Simulate event
    if idx_evt<=5       % cell division

        % Draw cell
        phn_i = idx_evt; % identify phenotype of cell
        cell_i = randsample(n_curr,1,true,phn_curr==phn_i);
        pos_i = pos_curr(cell_i,:);

        % Define hollow 'cubic sphere' with radius k around cell i and scan for
        % free spaces (pos_f)
        pos_f = NaN;
        k = 1;
        prevscan = [0 0 0];

        while isnan(sum(pos_f)) && k<=bdgmax

            % Begin with cube of length 2k+1
            nghbrd_k = combvec(-k:k,-k:k,-k:k)';

            % Remove previously scanned positions and those more than k
            % units from position i
            nghbrd_k(ismember(nghbrd_k,prevscan,'rows'),:) = [];

            euc_k = sqrt(sum(nghbrd_k.^2,2));
            nghbrd_k(euc_k>k,:) = [];

            % Scan for unoccupied positions k units from position i
            nghbrs_ik = pos_i+nghbrd_k;
            free_k = nghbrs_ik(~ismember(nghbrs_ik,pos_curr,'rows'),:);
            free_ik = free_k - pos_i;

            if isempty(free_ik)
                prevscan = [prevscan; nghbrd_k];
                k = k+1;                    
            else
                % Choose free position f, favoring those nearest center
                % of spheroid
%                     pos_cntr = mean(pos_curr);
%                     euc_cntrk = sqrt(sum((free_k-pos_cntr).^2,2));
%                     pos_f = free_k(randsample(size(free_k,1),1,true,exp(-5*euc_cntrk)),:);

                % Choose free position f closest to dividing cell i
                norms_ik = sqrt(sum(free_ik.^2,2));
                min_free_k = find(norms_ik==min(norms_ik));
                if length(min_free_k)==1
                    pos_f = free_k(min_free_k,:);
                else % Given the option, choose position f to be candidate space closest to center
                    pos_cntr = mean(pos_curr);
                    nrst_k = free_k(min_free_k,:);
                    dst_cntr_nrstk = sqrt(sum((nrst_k-pos_cntr).^2,2));
                    [~,min_dcn_idx] = min(dst_cntr_nrstk);
                    pos_f = nrst_k(min_dcn_idx,:);
                end                    
            end

        end

        % If unable to find free position within bdgmax distance,
        % continue to next iteration
        if isnan(pos_f)
            continue
        end

        % Push cells between i and f until free space adjacent to i
        % Definitions:
        % position i = site of dividing cell
        % position f = site of initially free space
        % position g = site of cell that moves to position f
        % position j = site of new cell (adjacent to position i)

        nghbrsi1 = pos_i + nghbrs1;
        nghbrsi1(ismember(pos_i,nghbrsi1,'rows'),:) = [];

        while ~ismember(pos_f,nghbrsi1,'rows')

            % Find position g adjacent to f closest to line through i
            % and f

            vec_if = pos_f - pos_i;
            dstr = NaN(26,1);
            for r = 1:26
                % Candidate for migrant cell g
                pos_gr = pos_f+nghbrs1(r,:);
                vec_ig = pos_gr-pos_i;

                % Consider only positions g that are closer to i than f
                % is to i
                if norm(vec_ig) < norm(vec_if)
                    % Use parallelogram definition of cross product to
                    % find distance between candidate position g and
                    % line passing through i and f
                    dstr(r) = norm(cross(vec_ig,vec_if))/norm(vec_if);
                end
            end

            % Choose position g
            idx_mindst = find(dstr==min(dstr));
            if length(idx_mindst)>1
                pos_g = pos_f+nghbrs1(randsample(idx_mindst,1),:);
            else
                pos_g = pos_f+nghbrs1(idx_mindst,:);
            end

            % Move cell in position g to position f
            pos_curr(ismember(pos_curr,pos_g,'rows'),:) = pos_f;

            % Position g is now open to receive next pushed cell
            pos_f = pos_g;

        end

        % Once a free space is adjacent to cell i, it divides ==> cell j
        pos_j = pos_f;
        pos_curr = [pos_curr; pos_j];

        % Identify phenotype of cells after budging
        phn_curr = [phn_curr; phn_i]; % assign same phenotype to daughter cells

        j = j+1;
        pos_next = pos_curr;
        phn_next = phn_curr;

%             disp(length(phn_curr));

    elseif idx_evt>5 && idx_evt<=25   % cell death

        % Draw cell
        phn_i = rem(idx_evt,5); % phenotype of cell is equal to remainder of idx_evt
        if phn_i==0
            phn_i = 5;
        end

        cell_i = randsample(n_curr,1,true,phn_curr==phn_i);


        % If dying cell_i was a gene drive cell, subtract one from counter
        % in appropriate rows of bystndr_mat
        if phn_i==3 && swtch2
            pos_i = pos_curr(cell_i,:);
            adj_i = pos_i+kill_rad_mat; % identify sites within kill radius of dying gene drive cell

            % Update bystndr_mat
            bystndr_adj_idx = ismember(bystndr_mat(:,1:3),adj_i,'rows');
            bystndr_mat(bystndr_adj_idx,4) = bystndr_mat(bystndr_adj_idx,4) - 1; % subtract 1 from count of associated gene drive cells for appropriate sites
            bystndr_mat(bystndr_mat(:,4)<0,4) = 0; % because this approach is not exact (bystndr_mat is calculated exactly ever bystndr_intrvl iterations, not every iteration), it's conceivable that the counter drops below zero; this fixes any issues that may cause
        end

        % Update state
        pos_curr(cell_i,:) = [];
        phn_curr(cell_i) = [];

        j = j+1;
        pos_next = pos_curr;
        phn_next = phn_curr;

%             disp(length(phn_curr));

    else % bystander effect clearance event

        zeros_state_idx = find(bystndr_mat(:,4)==0);

        row_clr = randsample(zeros_state_idx,1); % randomly select residual bystander effect row to clear
        bystndr_mat(row_clr,:) = []; % remove row

    end
    
    % Store state every str_intrvl
    if rem(j,str_intrvl)==1
        s = s+1;

        output(s,:) = [j t_curr pop_curr];
        
%         % Update complete state
%         state_store = NaN(nmax,4);
%         state_store(1:length(phn_curr),:) = [pos_curr phn_curr];
%         state(:,:,s) = state_store;
       
        disp([size0 t_curr n_curr])
    end

    % Next state becomes "current" state for next iteration
    pos_curr = pos_next;
    phn_curr = phn_next;
    n_curr = length(phn_curr);
    t_curr = sum(t);

    % Check if gene drive population has been depleted and, if so, if 10
    % days have passed
    if sum(phn_curr==3)==0 && ~gddepleted
        gddepleted = true;
        t_gddepleted = t_curr;
    end

    if gddepleted
        t_since_gddepleted = t_curr - t_gddepleted;
        if t_since_gddepleted >= 150
            gddepleted150 = true;
        end
    end

end

% Save final state

s = s+1;

% Current populations
nS = sum(phn_curr==1);
nR = sum(phn_curr==2);
nG = sum(phn_curr==3);
nAS = sum(phn_curr==4);
nAR = sum(phn_curr==5);
pop_curr = [nS nR nG nAS nAR];

output(s,:) = [j t_curr pop_curr];

% % Update complete state
% state_store = NaN(nmax,4);
% state_store(1:length(phn_curr),:) = [pos_curr phn_curr];
% state(:,:,s) = state_store;

      
% save('ExampleSimulation_size1e4_bdgmax30_091220.mat')
% 
% if sum(phn_curr==2)==0
%     t_100 = Inf;
% elseif size(state_curr,1) >= nmax
%     t_100 = NaN;
% else
%     t_100 = sum(t);
% end
% 
% %     disp(i)
% % output(i,:) = [bdgmax ncaf0 i t_100];
% 
% time = cumsum([0 t]);
% 
% phn = state(:,4,:);
% St = reshape(sum(phn==1),[j 1]);
% Rt = reshape(sum(phn==2),[j 1]);
% At = reshape(sum(phn==3),[j 1]);

% output = [time' St Rt At];
    
end


%% Plotting
% 
% % Plot population dynamics
% 
% figure
% plot(output(:,2),output(:,3),'Color',clr_sen)
% hold on
% plot(output(:,2),output(:,4),'Color',clr_res)
% plot(output(:,2),output(:,5),'Color',clr_gdv)
% plot(output(:,2),output(:,6),'Color',clr_ads)
% plot(output(:,2),output(:,7),'Color',clr_adr)
% hold off
% 
% 
% % Plot 3D tumor
% 
% minx = min(min(state(:,1,:)));
% maxx = max(max(state(:,1,:)));
% miny = min(min(state(:,2,:)));
% maxy = max(max(state(:,2,:)));
% minz = min(min(state(:,3,:)));
% maxz = max(max(state(:,3,:)));
% maxdim = max(abs([minx maxx miny maxy minz maxz]))+1;
% 
% figure('Position',[10 10 600 600])
% for n = 1:100
% for i = 1:size(state,3)
%     
%     Srows = state(:,4,i)==1;
%     Rrows = state(:,4,i)==2;
%     Grows = state(:,4,i)==3;
%     ASrows = state(:,4,i)==4;
%     ARrows = state(:,4,i)==5;
%     
%     clrs = NaN(nmax,3);
%     clrs(Srows,:) = repmat(clr_sen,sum(Srows),1);
%     clrs(Rrows,:) = repmat(clr_res,sum(Rrows),1);
%     clrs(Grows,:) = repmat(clr_gdv,sum(Grows),1);
%     clrs(ASrows,:) = repmat(clr_ads,sum(ASrows),1);
%     clrs(ARrows,:) = repmat(clr_adr,sum(ARrows),1);
% 
%     scatter3(state(:,1,i),state(:,2,i),state(:,3,i),3e3,clrs,'.');
%     hold on
%     scatter3(state(:,1,i),state(:,2,i),repmat(-maxdim,size(state,1),1),3e3,[.9 .9 .9],'.'); % add shadow
%     xlim([-maxdim maxdim]);
%     ylim([-maxdim maxdim]);
%     zlim([-maxdim maxdim]);
%     view(-75+i/5,20)
%     if output(i,2) < t_swtch2
%         title(strcat('Time: ',num2str(round(output(i,2),0)),' Days'),'Switch 1','Color',[249 154 0]/255)
%     else
%         title(strcat('Time: ',num2str(round(output(i,2),0)),' Days'),'Switch 2','Color',[186 80 255]/255)
%     end
%     grid off
%     set(gca,'TickLength',[0 0],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);%,'visible','off')
%     hold off
%     pause(1e-2)
%     
%     F(i) = getframe(gcf);
%     drawnow    
% end
% end
% 
% % create the video writer with 1 fps
% writerObj = VideoWriter('ExampleSim_091120b_bdgmax30.mp4','MPEG-4');
% writerObj.FrameRate = 20;
% 
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

% time = cumsum([0 t]);
% 
% phn = state(:,4,:);
% St = reshape(sum(phn==1),[j 1]);
% Rt = reshape(sum(phn==2),[j 1]);
% At = reshape(sum(phn==3),[j 1]);
% 
% figure
% plot(time,[St Rt At])
% xlabel('Days')
% ylabel('Cell Number')
% legend({'sen' 'res' 'adj'},'Location','northeastoutside')
