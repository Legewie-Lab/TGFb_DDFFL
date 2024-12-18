function setup_smadffl(vl)
% Function to set up and run SMAD feed-forward loop analysis

% Configure parallel cluster settings
pc = parcluster('local');

% Get the number of dedicated cores from the SLURM environment
num_workers = str2num(getenv('SLURM_NPROCS'));

% Set up a temporary directory for parallel processing
parpool_tmpdir = [getenv('TMP'), '/.matlab/local_cluster_jobs/slurm_jobID_', getenv('SLURM_JOB_ID')];
mkdir(parpool_tmpdir);
pc.JobStorageLocation = parpool_tmpdir;

% Start the parallel pool
parpool(pc, num_workers);

% Load the gene names and dataset
load gename.mat; % Gene names
load('alles_d2d.mat'); % Dataset with statistical analysis results

% Filter data with p-value < 0.05
rej_da = alles_d2d(round(alles_d2d{:,'pvalue'}, 2) < 0.05, :);

% Extract gene IDs
nam = rej_da{:,'geneID'};

% Add the required Data2Dynamics library to the MATLAB path
addpath(genpath('/home/st/st_us-040512/st_ac137903/smad/D2D/Data2Dynamics-d2d-28a28b1'));
savepath pathsa.m;

% Define the layout for the analysis grid
numRows = 49; % Number of rows
numCols = 36; % Number of columns

% Create a sequential vector and reshape it into a grid
vec = 1:(numRows * numCols);
A_ma = reshape(vec, numCols, numRows)';

% Specify model and number of fitting iterations
model = [1];
fit = 500; % Latin Hypercube Sampling iterations


% Run the analysis for each gene in parallel
parfor i = 1:numel(nam)
    smadfitffl(nam{i}, model, fit);
end

end
