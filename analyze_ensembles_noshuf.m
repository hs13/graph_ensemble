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

% %% crf_Rz1pt5_noshuf_
% mousedate = 'MU21_1/220113/';
% mousedate = 'MU21_2/220110/'; % 0 IC-encoders and 2 RC-encoders in this session -- maybe skip
% mousedate = 'HS_CamKIIGC6s_51/210225/';
% mousedate = 'MU21_2/220113/';

mdsplit = strsplit(mousedate, '/');
mdcond = [mdsplit{1} '_' mdsplit{2}];
datadir = '/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed/graph_ensemble_data/';
noshufdir = ['/Users/hyeyoung/Documents/DATA/graph_ensemble_results/Rz1pt5_noshuf_' mdcond '/'];
pathpp = ['/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed/' mousedate];
disp(mousedate)

load([datadir 'Rz1pt5_noshuf_' mdcond '.mat']); % Loads variables `data` and `stimuli`
best_model = load([noshufdir 'results/best_model_full.mat']);
% shuffle_model = load([resultdir 'shuffled_experiment_demo/results/fulldata.mat']);

%% find_temporal_ens_nodes_noshuf is called by find_plot_temporal_crf_core_noshuf
if exist([noshufdir 'crf_Rz1pt5_noshuf_' mdcond '.mat'], 'file')
load([noshufdir 'crf_Rz1pt5_noshuf_' mdcond '.mat'])
else
% [ens_nodes, results0] = find_temporal_ens_nodes_noshuf(best_model, data, stimuli);

% took 42min for MU21_1_220113 (572 neurons, 2800 stimuli)
tic
[ens_nodes, crfresults] = find_plot_temporal_crf_core_noshuf(best_model,data,stimuli);
toc

save([pathpp 'crf_Rz1pt5_noshuf.mat'], 'ens_nodes', 'crfresults')
save([noshufdir 'crf_Rz1pt5_noshuf_' mdcond '.mat'], 'ens_nodes', 'crfresults')
end

%%
[num_frame, num_stim] = size(stimuli);
num_node = size(best_model.graph,1);
num_orig_neuron = size(data, 2);
time_span = best_model.time_span;

core_crf = crfresults.core_crf;
epsum = crfresults.epsum;
auc = crfresults.auc;
auc_ens = crfresults.auc_ens;

neurongroups = load(strcat(pathpp, 'neurongroups.mat'));
load(strcat(pathpp, 'ICsig.mat'))
ICtt = [0 106 107 110 111 1105 1109];

%% plot
nodesz = 30;
nsmi = min(epsum);
nsma = max(epsum);
aucmi = 0.25;%0;
aucma = 0.75;%1;

f = figure; set(gcf,'color','w')
f.Name = sprintf('K=%d', time_span);
color_by_offset = @(x) floor((x-1)/num_orig_neuron) / max(1, time_span-1);
for ii = 2:5%1:num_stim
    
    % AUC - node strength plot
    %cur_axes = subplot(Nrows, Ncols,ii); hold on
    cur_axes = subplot(2, 2,ii-1); hold on
    
    colormap(cur_axes, autumn)
    scatter(epsum,auc(:,ii),nodesz,0.5*[1 1 1],'filled')
    % Stimuli nodes blue
    scatter(epsum(end - num_stim + 1:end),auc(end - num_stim + 1:end,ii),nodesz,[0 0 1],'filled')
    % Core nodes colored red->yellow according to how frame-offset
    scatter(epsum(core_crf{ii}),auc(core_crf{ii},ii),nodesz,arrayfun(color_by_offset, core_crf{ii}),'filled')
    % Active stimulus node cyan
    scatter(epsum(num_node - num_stim + ii),auc(num_node - num_stim + ii,ii),nodesz,[0 1 1],'filled')

    % IC-encoders/RC-encoders    
    ttcol = [0 0.5 0; 0.8 0.4 0; 1 0.6 0.2; 0.2 1 0.2];
    for isp =1:4
        switch isp
            case 1
                neuoi = 'ICencoder1';
            case 2
                neuoi = 'RCencoder1';
            case 3
                neuoi = 'RCencoder2';
            case 4
                neuoi = 'ICencoder2';
        end
        scatter(epsum(neurongroups.(neuoi)),auc(neurongroups.(neuoi),ii),nodesz,'filled','MarkerFaceColor', ttcol(isp,:))
    end
    scatter(epsum(inducerencoder),auc(inducerencoder,ii),nodesz,'filled','MarkerFaceColor', [0 0 0])
    
    plot([nsmi nsma],mean(auc_ens{ii})*[1 1],'k--');
    plot([nsmi nsma],(mean(auc_ens{ii})+std(auc_ens{ii}))*[1 1],'--',...
        'color',0.7*[1 1 1]);
    plot([nsmi nsma],(mean(auc_ens{ii})-std(auc_ens{ii}))*[1 1],'--',...
        'color',0.7*[1 1 1]);
    plot(nanmean(epsum)*[1 1],[aucmi aucma],'k--');
%     plot((nanmean(epsum)+nanstd(epsum)/10)*[1 1],[aucmi aucma],'--',...
%         'color',0.7*[1 1 1]);
%     plot((nanmean(epsum)-nanstd(epsum)/10)*[1 1],[aucmi aucma],'--',...
%         'color',0.7*[1 1 1]);
    xlim([nsmi nsma]); ylim([aucmi aucma])
    xlabel('node strength'); ylabel(['AUC' num2str(ICtt(ii))]);
    title(['Trial Type ' num2str(ICtt(ii))], 'Color', ttcol(ii-1,:))
        
end

%%
Nneurons = length(inducerencoder);
encgcols = [.5 .5 .5; 0 0.7 0; 1 0.5 0; 0 0 0];
Nneugcum = [];
neugcum = [];
epsumgcum = [];
aucgcum = [];
ytl = {};
for isp =1:4
    switch isp
        case 1
            neuoi = true(Nneurons, 1);
            neutit = 'All';
        case 2
            neuoi = ICencoder;
            neutit = 'IC-enc';
        case 3
            neuoi = RCencoder;
            neutit = 'RC-enc';
        case 4
            neuoi = inducerencoder;
            neutit = 'ind-enc';
    end
    if nnz(neuoi)==0
        continue
    end
    ytl = cat(1,ytl,neutit);
    Nneugcum = cat(1, Nneugcum, nnz(neuoi));
    neugcum = cat(1, neugcum, isp*ones(nnz(neuoi),1) );
    epsumgcum = cat(1, epsumgcum, epsum(neuoi));
    aucgcum = cat(1, aucgcum, auc(neuoi,:));
end

[P,ANOVATAB,STATS] = kruskalwallis(epsumgcum, neugcum);
[c,m,h] = multcompare(STATS);
set(gca, 'YTick', h.CurrentAxes.YTick, 'YTickLabel', flip(ytl), 'FontSize', 11)
xlabel('Node Strength');
set(h, 'Position',[300 200 400 100])

fs=16;
figure('Position',[300*(isp-1), 100, 300 330])
t = tiledlayout(1,1);
ax1 = axes(t);
xl = [0.5 length(ytl)+0.5];
hold on
b = boxplot(epsumgcum, neugcum, 'Colors', encgcols, 'notch', 'on', 'Symbol', 'k+');
set(b,{'linew'},{3})
yl = ylim;
text(xl(2), yl(1)+0.01*range(yl), sprintf('Pkw=%.4f', P), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs)
xlim(xl)
% if doyl; ylim(yl); end
set(ax1, 'FontSize', fs, 'XTick', 1:length(ytl),  'XTickLabel', ytl, 'XTickLabelRotation', 45)
% xlabel(ax1, 'Neuron Subset')
ylabel(ax1, 'CRF Node Strength')

ax2 = axes(t);
% ax2.Xlim = xl;
xlim(ax2, xl)
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.YAxis.Visible = 'off';
ax2.Color = 'none';
xtl2 = strsplit(sprintf('%d\n', Nneugcum), '\n');
set(ax2, 'FontSize', fs, 'XTick', 1:length(ytl),  'XTickLabel', xtl2(1:end-1) )
xlabel(ax2, '# Neurons in Subset')

%%
Nneurons = length(inducerencoder);
enccols = [.5 .5 .5; 0 0.5 0; 0.8 0.4 0; 1 0.6 0.2; 0.2 1 0.2; 0 0 0];
Nneucum = [];
neucum = [];
epsumcum = [];
auccum = [];
ytl = {};
for isp =1:6
    switch isp
        case 1
            neuoi = true(Nneurons, 1);
            neutit = 'All';
        case 2
            neuoi = neurongroups.ICencoder1;
            neutit = 'IC1-enc';
        case 3
            neuoi = neurongroups.RCencoder1;
            neutit = 'RC1-enc';
        case 4
            neuoi = neurongroups.RCencoder2;
            neutit = 'RC2-enc';
        case 5
            neuoi = neurongroups.ICencoder2;
            neutit = 'IC2-enc';
        case 6
            neuoi = inducerencoder;
            neutit = 'ind-enc';
    end
    if nnz(neuoi)==0
        continue
    end
    ytl = cat(1,ytl,neutit);
    Nneucum = cat(1, Nneucum, nnz(neuoi));
    neucum = cat(1, neucum, isp*ones(nnz(neuoi),1) );
    epsumcum = cat(1, epsumcum, epsum(neuoi));
    auccum = cat(1, auccum, auc(neuoi,:));
end

[P,ANOVATAB,STATS] = kruskalwallis(epsumcum, neucum);
[c,m,h] = multcompare(STATS);
set(gca, 'YTick', h.CurrentAxes.YTick, 'YTickLabel', flip(ytl), 'FontSize', 11)
xlabel('Node Strength');
set(h, 'Position',[300 200 400 100])

% ii=5;
% [P,ANOVATAB,STATS] = kruskalwallis(auccum(:,ii), neucum);
% [c,m,h] = multcompare(STATS);
% set(gca, 'YTick', h.CurrentAxes.YTick, 'YTickLabel', flip(ytl), 'FontSize', 11)
% xlabel(['AUC' num2str(ICtt(ii))]);

%%
fs=16;
figure('Position',[300*(isp-1), 100, 300 330])
t = tiledlayout(1,1);
ax1 = axes(t);
xl = [0.5 6.5];
hold on
b = boxplot(epsumcum, neucum, 'Colors', enccols, 'notch', 'on', 'Symbol', 'k+');
set(b,{'linew'},{3})
yl = ylim;
text(xl(2), yl(1)+0.01*range(yl), sprintf('Pkw=%.4f', P), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs)
xlim(xl)
% if doyl; ylim(yl); end
set(ax1, 'FontSize', fs, 'XTick', 1:length(ytl),  'XTickLabel', ytl, 'XTickLabelRotation', 45)
% xlabel(ax1, 'Neuron Subset')
ylabel(ax1, 'CRF Node Strength')

ax2 = axes(t);
% ax2.Xlim = xl;
xlim(ax2, xl)
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.YAxis.Visible = 'off';
ax2.Color = 'none';
xtl2 = strsplit(sprintf('%d\n', Nneucum), '\n');
set(ax2, 'FontSize', fs, 'XTick', 1:length(ytl),  'XTickLabel', xtl2(1:end-1) )
xlabel(ax2, '# Neurons in Subset')

%%
figure('Position', [100 200 1200 300])
for ii = 2:5%1:num_stim
    subplot(1, 4, ii-1)
xl = [0.5 6.5];
hold on
b = boxplot(auccum(:,ii), neucum, 'Colors', enccols, 'notch', 'on', 'Symbol', 'k+');
set(b,{'linew'},{3})
yl = ylim;
[P,ANOVATAB,STATS] = kruskalwallis(auccum(:,ii), neucum, 'off');
text(xl(2), yl(1)+0.01*range(yl), sprintf('Pkw=%.2f', P), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs)
xlim(xl)
% if doyl; ylim(yl); end
set(gca, 'FontSize', fs, 'XTick', 1:length(ytl),  'XTickLabel', ytl, 'XTickLabelRotation', 45)
% xlabel(ax1, 'Neuron Subset')
ylabel(['CRF AUC' num2str(ICtt(ii))]);
end
