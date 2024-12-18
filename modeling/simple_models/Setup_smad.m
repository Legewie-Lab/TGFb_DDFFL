% Load the list of gene names
load gename.mat; % 'gename' contains the list of gene names to analyze

% Specify the model type to use
model = 2; % Set to 2 for inhibitor model
% model = 1; % Uncomment to use activator model

% Define the number of Latin Hypercube Sampling (LHS) fits
fit = 500; % Number of LHS iterations for parameter fitting

% Define the indices of genes to analyze
data_index = [1:length(gename)]; % Analyze all genes in the list
% data_index = [2396 4807 720 4104]; % Uncomment and modify to analyze specific genes by their indices

% Parallel loop to fit models for each gene
parfor i = data_index
    smadfit(gename{i}, model, fit); % Perform model fitting for the current gene
end

% Uncomment the line below to run the `smadfit` function for a single gene
% smadfit(name, model, fit);