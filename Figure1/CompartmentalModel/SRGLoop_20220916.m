%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop through SRG simulations %%%
%%%    Scott Leighow 09/22/21    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close

rng(10)

% parpool(48)

% might change default mut to 1e-8 and pop_treat to 1e9 (a bit more
% reasonable, and then sensitivity analysis spans more on both sides of
% default parameters)

%% Parameters

npats = 48;   % number of patients to be simulated per parameter set

pop_treat_vec = 1e9; %logspace(8,12,5); % Default 1e9
mut_vec = 1e-8;      % Default 1e-8 
inf_eff_vec = [0 logspace(-3,-0.5,6)];
swtch_dly_vec = 1*365;
dim_res_vec = 0:0.02:0.14;
a_B_vec = 0.04;
trnovr_vec = 14; %[1 7 14 21]; % Default 14
ngr_vec = [0.05 0.01 0.015 0.02 0.025]; % Default 0.01

%% Loop Through Simulations

for pop_treat = pop_treat_vec
    for mut = mut_vec
        for inf_eff = inf_eff_vec
            for swtch_dly = swtch_dly_vec
                for dim_res = dim_res_vec
                    for a_B = a_B_vec
                        for trnovr = trnovr_vec
                            for ngr = ngr_vec
                                
                                b = ngr*trnovr;
                                d = ngr*(trnovr-1);
            
                                out = SRG_20220916(pop_treat,mut,inf_eff,swtch_dly,dim_res,a_B,b,d,npats);

                                fname = strcat('SRG_logpoptreat',num2str(log10(pop_treat)), ...
                                              '_neglogmut',num2str(-log10(mut)), ...
                                              '_negloginfeff',num2str(-log10(inf_eff)), ...
                                              '_swtchdly',num2str(swtch_dly), ...
                                              '_dimres',num2str(dim_res), ...
                                              '_aB',num2str(a_B), ...
                                              '_turnovr',num2str(trnovr), ...
                                              '_ngr',num2str(ngr), ...
                                              '_20220916.csv');

                                csvwrite(fname,out);

                                disp([ngr -log10(inf_eff) dim_res])
                                
                            end
                        end
                    end
                end
            end
        end
    end
end