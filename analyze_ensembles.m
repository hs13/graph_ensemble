%{
Finding core ensembles in the model
With a trained CRF model on the dataset of interest, and a collection of 
models trained on shuffled versions of the dataset, use 
src/core/find_temporal_ens_nodes.m to find the ensembles corresponding to 
each stimulus. We will continue our previous example:

Start matlab. You can do so from the terminal with the following command:
matlab -nodesktop -nosplash -nodisplay
Load the models and data in matlab:
addpath(genpath(‘~/graph_ensemble/’)); % This should be your source_directory
cd expt/
best_model = load('experiment_demo/results/best_model_full.mat');
shuffle_model = load('shuffled_experiment_demo/results/fulldata.mat');
load('~/data/experiment_demo.mat'); % Loads variables `data` and `stimuli`

Find ensemble nodes:
ens_nodes = find_temporal_ens_nodes(best_model, shuffle_model, data, stimuli)

ens_nodes is a cell vector where each cell contains the ensemble neurons 
found for each stimuli. Each such stimuli cell contains a further cell 
vector where each cell contains the ensemble neurons found for each offset 
frame of the time_span window, with the first corresponding to no offset.

Another script, scripts/core/find_plot_temporal_crf_core.m, can also be 
used on a desktop system to find ensemble neurons and plot some features, 
including spatial arrangement if coordinates are provided.
%}
addpath(genpath('/Users/hyeyoung/Documents/CODE/graph_ensemble/'))

mdcond = 'MU21_1_220113';
datadir = '/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed/graph_ensemble_data/';
resultdir = ['/Users/hyeyoung/Documents/DATA/graph_ensemble_results/Rz1pt5_' mdcond '/'];
shufdir = ['/Users/hyeyoung/Documents/DATA/graph_ensemble_results/shuffled_Rz1pt5_' mdcond '/'];

% best_model = load('experiment_demo/results/best_model_full.mat');
% shuffle_model = load('shuffled_experiment_demo/results/fulldata.mat');
% load([sourcedir 'experiment_demo.mat']); % Loads variables `data` and `stimuli`

load([datadir 'Rz1pt5_noshuf_' mdcond '.mat']); % Loads variables `data` and `stimuli`
best_model = load([resultdir 'results/best_model_full.mat']);
shuffle_model = load([shufdir 'results/fulldata.mat']);

% the following is called by find_plot_temporal_crf_core
% [ens_nodes, results0] = find_temporal_ens_nodes(best_model, shuffle_model, data, stimuli);

tic
[ens_nodes, results] = find_plot_temporal_crf_core(best_model,shuffle_model,data,stimuli, coords);
toc

% % what perceptage of frames is each neuron active
% disp(mean(data,1))
% mean(mean(data,1)) % 0.096
% median(mean(data,1)) % 0.07

