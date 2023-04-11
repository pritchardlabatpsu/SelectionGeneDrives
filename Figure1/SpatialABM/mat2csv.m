%%% Rewrite all mat files in folder to csv files

clear

%% Get files
fileinfo = dir('*.mat');
fnames = {fileinfo.name};

%% Write mat files as csv

for i = 1:length(fnames)

    fname_i = fnames{i};    % select ith file
    load(fname_i)           % read in data data

    fname_csv = regexprep(fname_i,"mat","csv");
    csvwrite(fname_csv,output); % save file

    clear output % clear output variable to make sure it's called fresh each iteration
    delete(fname_i) % delet .mat file

end
