
% Load the data table
% Note: The file path needs to point to the correct directory.
tabel=readtable('../../4823_DEG_upon_low_high_TGFb_dose.csv');

% Load the SMAD signal data
smadsig = readtable('data_smad.xls');

% Remove duplicate rows from the data table
tabel = unique(tabel);

% Define column indices for low and high values
low = [3:8];
high = [9:14];

% Extract low and high data as transposed matrices
low_data = tabel{:, low}';
high_data = tabel{:, high}';

% Define the time steps for analysis and convert to a table
time_step = [45; 90; 180; 360; 720; 1440];
time_tab = array2table(time_step, 'VariableNames', {'time'});

% Clean up gene names by removing special characters
gename = tabel{:, end};
for i = 1:numel(gename)
    gename{i} = strrep(gename{i}, '-', '');    
    gename{i} = strrep(gename{i}, '.', '');
end

% Convert SMAD signal data to an array
smadsig_array = table2array(smadsig);

% Preallocate arrays for modified high and low data
high_mod = zeros(length(smadsig_array) + 1, length(high_data));
low_mod = zeros(length(smadsig_array) + 1, length(high_data));

% Assign the last row of the original high and low data
high_mod(end, :) = high_data(end, :);
low_mod(end, :) = low_data(end, :);

% Find indices of time steps in the SMAD signal array
[~, time_idx] = ismember(time_step, smadsig_array(:, 1));
time_idx(end, :) = [];

% Populate high and low mod arrays based on time indices
high_mod(time_idx, :) = high_data(1:5, :);
low_mod(time_idx, :) = low_data(1:5, :);

% Replace zero values with empty entries
high_mod = cellfun(@(x) x(logical(x)), num2cell(high_mod), 'uni', false);
low_mod = cellfun(@(x) x(logical(x)), num2cell(low_mod), 'uni', false);

% Generate gene names for high and low data
gename_high = strcat(gename, '_high');
gename_low = strcat(gename, '_low');

% Convert high and low mod arrays to tables
high_mod = array2table(high_mod, 'VariableNames', gename_high);
low_mod = array2table(low_mod, 'VariableNames', gename_low);

% Update SMAD signal array and table to match the data structure
smadsig_array(length(smadsig_array) + 1, :) = [1440, 0];
smadsig = cell2table(cellfun(@(x) x(logical(x)), num2cell(smadsig_array), 'uni', false), 'VariableNames', {'time', 'SMADx'});
smadsig.time{1} = 0;
smadsig.SMADx{1} = 0;

% Create a table for high data
data_table_high = array2table(high_data, 'VariableNames', gename);

% Write CSV files for each gene
parfor i = 1:size(data_table_high, 2)
    % Create a combined table with SMAD signal and gene expression data
    T = table(smadsig_array(:, 1), high_mod{:, i}, low_mod{:, i});
    T.Properties.VariableNames = ["time", "gnexhigh_au", "gnexlow_au"];
    writetable(T, ['Data/data_' gename{i} '.csv']);
end

% Copy and modify template definition files for each gene
parfor i = 1:size(data_table_high, 2)
    copyfile('template_new.def', ['Data/data_' gename{i} '.def']);
    fileID = fopen(['Data/data_' gename{i} '.def'], 'a');
    fprintf(fileID, '\n sd_gnexhigh    "0.11*gnexhigh_au"');
    fprintf(fileID, '\n sd_gnexlow    "0.11*gnexlow_au"');
    fclose(fileID);
end

