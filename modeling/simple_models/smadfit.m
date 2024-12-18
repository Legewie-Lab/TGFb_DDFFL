function smadfit(name, model, fit)
% Function to fit SMAD signaling model parameters to data.
% Inputs:
%   - name: Gene name to analyze.
%   - model: Model type (1 for activator, 2 for inhibitor).
%   - fit: Number of LHS fits to perform.

% Close all figures and initialize the environment
close all
arInit;

% Load gene names
load gename.mat; % gename contains a list of gene names

% Find the index of the gene of interest
ind = find(ismember(gename, name));

% Load the appropriate model based on the input parameter
if model == 1
    arLoadModel('Simple_mode_SMAD_activ'); % Load activator model
elseif model == 2
    arLoadModel('Simple_mode_SMAD_inhib'); % Load inhibitor model
end

% Load data for the selected gene
arLoadData(['data_' gename{ind}]);

% Compile the model and data
arCompileAll;

% Set initial parameters with bounds and log10 scale
arSetPars('k', 1, 1, 1, -3, 3);  % Parameter k
arSetPars('l', -1, 1, 1, -5, 1); % Parameter l
arSetPars('r', -3, 1, 1, -5, 2); % Parameter r
arSetPars('h', 1, 1, 0, 0, 10);  % Parameter h (fixed to 1)

% Set the standard deviation for gene expression data fitting
arSetPars('sd_gnex', 0.11, 0, 0, -5, 5);

% Configure error model (no error model in this case)
ar.config.fiterrors = 0;

% Perform model fitting
arFit;

% Perform Latin Hypercube Sampling (LHS) fitting
arFitLHS(fit);

% Extract fitted parameters and statistics
dr = ar.model.data;  % Access model data structure
ksolu = dr.pNum;     % Extract parameter values
kname = dr.p;        % Extract parameter names

% Define time points and corresponding data labels
time_step = [{'45_high'}, {'90_high'}, {'180_high'}, {'360_high'}, {'720_high'}, {'1440_high'}, ...
             {'45_low'}, {'90_low'}, {'180_low'}, {'360_low'}, {'720_low'}, {'1440_low'}];

time_data = [{'45_highdata'}, {'90_highdata'}, {'180_highdata'}, {'360_highdata'}, {'720_highdata'}, {'1440_highdata'}, ...
             {'45_lowdata'}, {'90_lowdata'}, {'180_lowdata'}, {'360_lowdata'}, {'720_lowdata'}, {'1440_lowdata'}];

% Create table header for storing results
tab_title = [time_step, kname, 'chisquare_high', 'chisquare_low', 'pvalue', time_data];

prd = []; % Predicted data
dta = []; % Observed data
chis = []; % Chi-squared values

% Extract predictions, observations, and chi-squared values
for i = 2:7
    prd(i-1, :) = [dr(i).tExp, dr(i).yExpSimu]; % Predicted data
    dta(i-1, :) = [dr(i).tExp, dr(i).yExp];     % Observed data
    chis(i-1, :) = dr(i).chi2;                  % Chi-squared values
end

chis = sum(chis, 1); % Sum chi-squared values

% Sort data by time points
dta = sortrows(dta, 1);
prd = sortrows(prd, 1);

% Extract values for data and predictions
val_data = (dta(:, 2:3))';
pred_activ = (prd(:, 2:3))';

% Generate filename based on model type
if model == 1
    string_name = ['fit_activ_' gename{ind}];
elseif model == 2
    string_name = ['fit_inhub_' gename{ind}];
end

% Perform chi-squared test and get p-value
arChi2Test;
pvl = ar.pval;

% Determine model classification based on p-value
if model == 1
    if pvl >= 0.05
        estm = {'activator'};
    else
        estm = {'rejected'};
    end
elseif model == 2
    if pvl >= 0.05
        estm = {'inhibitor'};
    else
        estm = {'rejected'};
    end
end

% Create a table of results
valus = array2table([pred_activ(1, :), pred_activ(2, :), ksolu, chis, pvl, val_data(1, :), val_data(2, :)], 'VariableNames', tab_title);

geneID = gename(ind);
genenm = table(geneID);

predic = table(estm); % Add prediction to the table

valus = [genenm, valus, predic];

eval([string_name ' = valus;']); % Save the table to a variable

% Save results and model states
if model == 1
    save(['model_fitting/fit_activ_' gename{ind} '.mat'], ['fit_activ_' gename{ind}]);
    arSave(['model_activ_' gename{ind} 'lhs_' num2str(fit) '.mat']);
elseif model == 2
    save(['model_fitting/fit_inhib_' gename{ind} '.mat'], ['fit_inhub_' gename{ind}]);
    arSave(['model_inhib_' gename{ind} 'lhs_' num2str(fit) '.mat']);
end

end
