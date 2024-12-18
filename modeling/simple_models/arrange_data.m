% Import the data tables
tabel=readtable('../../4823_DEG_upon_low_high_TGFb_dose.csv');
smadsig = readtable('data_smad.xls'); % SMAD signaling data

% Select column ranges for low and high expression values
low = [3:8]; % Columns corresponding to 'low' expression data
high = [9:14]; % Columns corresponding to 'high' expression data

% Extract data for the low and high value columns
low_data = tabel{:, low}';
high_data = tabel{:, high}';

% Define time steps for the data and convert to a table
time_step = [45; 90; 180; 360; 720; 1440]; % Time points in minutes
time_tab = array2table(time_step, 'VariableNames', {'time'}); % Time steps table

% Clean up gene names by removing '-' and '.'
gename = tabel{:, 3}; % Extract gene names
for i = 1:numel(gename)
    gename{i} = strrep(gename{i}, '-', '');    
    gename{i} = strrep(gename{i}, '.', '');
end

% Convert SMAD signaling data to an array for processing
smadsig_array = table2array(smadsig);

% Preallocate arrays for modified high and low data
high_mod = zeros(length(smadsig_array) + 1, length(high_data));
low_mod = zeros(length(smadsig_array) + 1, length(high_data));

% Fill the last row of modified data arrays with the corresponding high/low data
high_mod(end, :) = high_data(end, :);
low_mod(end, :) = low_data(end, :);

% Map time steps to indices in the SMAD signal array
[~, time_idx] = ismember(time_step, smadsig_array(:, 1)); % Find indices for time points
time_idx(end, :) = []; % Exclude the last entry (already handled)

% Assign data values for the matched time indices
high_mod(time_idx, :) = high_data(1:5, :);
low_mod(time_idx, :) = low_data(1:5, :);

% Replace zero values with empty cells
high_mod = cellfun(@(x) x(logical(x)), num2cell(high_mod), 'uni', false);  
low_mod = cellfun(@(x) x(logical(x)), num2cell(low_mod), 'uni', false);  

% Append '_high' and '_low' to gene names for high and low data tables
gename_high = strcat(gename, '_high');
gename_low = strcat(gename, '_low');

% Convert modified data arrays to tables
high_mod = array2table(high_mod, 'VariableNames', gename_high);
low_mod = array2table(low_mod, 'VariableNames', gename_low);

% Extend the SMAD signaling data to match other data structures
smadsig_array(length(smadsig_array) + 1, :) = [1440, 0]; % Add a time point at 1440 minutes
smadsig = cell2table(cellfun(@(x) x(logical(x)), num2cell(smadsig_array), 'uni', false), 'VariableNames', {'time', 'SMADx'});  

% Restore empty cells to zeros for consistency
smadsig.time{1} = 0;
smadsig.SMADx{1} = 0;

% Create a table for the high data with gene names as variable names
high_data_tabl = array2table(high_data, 'VariableNames', gename);

% Write data to CSV files for each gene using a parallel loop
parfor i = 1:size(high_data_tabl, 2)
    % Combine time, high, and low data into one table
    T = [smadsig(:, 1), high_mod(:, i), low_mod(:, i)];
    T.Properties.VariableNames = ["time", "gnexhigh_au", "gnexlow_au"]; % Set column names
    
    % Save the combined table to a CSV file named after the gene
    writetable(T, ['Data/data_' gename{i} '.csv']); 
end
