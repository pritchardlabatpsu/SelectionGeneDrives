%%% Call GeneDriveABMfunc

clear
close all
rng(100)

% Game plan:
% Previous results showed no change in outcome w/r/t initial seeding

% Parameters
size0_vec = 1e4;
tight_vec = [0 0.175 0.31 0.515 3]; % chosen to get roughly [1 0.75 0.5 0.25 0] of maximum dispersion for size0=1e4
n_foci_vec = 1;
kill_rad_vec = 1:5;
swtch_dly_vec = 180;
iter_vec = 6:20;

% Combine parameter vectors into combinatorial matrix
param_mat = combvec(size0_vec,tight_vec,n_foci_vec,kill_rad_vec,swtch_dly_vec,iter_vec)';

parfor n = 1:size(param_mat,1)

    % Parameters
    size0_n = param_mat(n,1);
    tight_n = param_mat(n,2);
    n_foci_n = param_mat(n,3);
    kill_rad_n = param_mat(n,4);
    swtch_dly_n = param_mat(n,5);
    iter_n = param_mat(n,6);

    % Call simulation function
    [output] = GeneDriveABMfunc_091422(size0_n,tight_n,n_foci_n,kill_rad_n,swtch_dly_n,iter_n);

    % Store results
    fname = strcat('GeneDriveABMResults_size0',num2str(size0_n), ...
                   '_tight',num2str(tight_n), ...
                   '_nfoci',num2str(n_foci_n),...
                   '_killrad',num2str(kill_rad_n),...
                   '_swtchdly',num2str(swtch_dly_n),...
                   '_iter',num2str(iter_n),...
                   '_091422.mat');

    parsave(fname,output)

end

function parsave(fname,output)
  save(fname, 'output')
end
