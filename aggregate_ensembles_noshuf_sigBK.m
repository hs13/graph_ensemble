addpath(genpath('/Users/hyeyoung/Documents/CODE/graph_ensemble/'))

% %% crf_noshuf_Rz1pt5_sigBK_
ses2agg = {'HS_CamKIIGC6s_51/210225/','HS_CamKIIGC6s_51/210226/','HS_CamKIIGC6s_53/210225/', ...
    'HS_CamKIIGC6s_53/210226/','HS_VIPHalo_4/210422/','HS_VIPHalo_4/210423/', ...
    'MU21_1/220110/','MU21_1/220113/','MU21_2/220110/','MU21_2/220111/','MU21_2/220113/'}; % 

neurongroupnames = {'indenc1', 'indenc2', 'indenc3', 'indenc4', ...
    'ICencoder1', 'ICencoder2', 'RCencoder1', 'RCencoder2', ...
    'indenc13', 'indenc24', 'indenc14', 'indenc23', ...
    'ICresp1', 'RCresp1', 'RCresp2', 'ICresp2', 'ICtuned1', 'RCtuned1', 'RCtuned2', 'ICtuned2'};
crfsigBK = struct();
crfsigBK.epsum = [];
crfsigBK.auc = [];
crfsigBK.zepsum = [];
crfsigBK.zauc = [];
crfsigBK.epmean = [];
crfsigBK.ICencoder = [];
crfsigBK.RCencoder = [];
crfsigBK.inducerencoder = [];
for g = 1:numel(neurongroupnames)
crfsigBK.(neurongroupnames{g}) = [];
end
Nneucum = zeros(length(ses2agg), 4);
epsumacc = NaN(length(ses2agg), 4);
aucacc = NaN(length(ses2agg), 4, 5);
zepsumacc = NaN(length(ses2agg), 4);
zaucacc = NaN(length(ses2agg), 4, 5);
epmeanacc = NaN(length(ses2agg), 4);
for ises = 1:numel(ses2agg)
mousedate = ses2agg{ises};

mdsplit = strsplit(mousedate, '/');
mdcond = [mdsplit{1} '_' mdsplit{2}];
datadir = '/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed/graph_ensemble_data/';
noshufdir = ['/Users/hyeyoung/Documents/DATA/graph_ensemble_results/noshuf_Rz1pt5_sigBK_' mdcond '/'];
pathpp = ['/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed/' mousedate];
disp(mousedate)

load([datadir 'noshuf_Rz1pt5_sigBK_' mdcond '.mat']); % Loads variables `data` and `stimuli`
best_model = load([noshufdir 'results/best_model_full.mat']);
% shuffle_model = load([resultdir 'shuffled_experiment_demo/results/fulldata.mat']);

%% find_temporal_ens_nodes_noshuf is called by find_plot_temporal_crf_core_noshuf
load([noshufdir 'crf_noshuf_Rz1pt5_sigBK_' mdcond '.mat']) % 'ens_nodes', 'crfresults'
neurongroups = load(strcat(pathpp, 'neurongroups.mat'));
load(strcat(pathpp, 'ICsig.mat'))
ICtt = [0 106 107 110 111 1105 1109];
neusigBK = find(PkwBK<0.05);

num_stim = size(crfresults.stimuli,2);
Nneurons = length(neusigBK);

epsum = crfresults.epsum(1:end-num_stim);
auc = crfresults.auc(1:end-num_stim,:);

crfsigBK.epsum = cat(1, crfsigBK.epsum, epsum );
crfsigBK.auc = cat(1, crfsigBK.auc, auc );

% epmean = epsum/Nneurons;
epmean = epsum/(Nneurons+num_stim);
crfsigBK.epmean = cat(1, crfsigBK.epmean, epmean );

zepsum = (crfresults.epsum(1:end-num_stim) - mean(crfresults.epsum(1:end-num_stim))) / std(crfresults.epsum(1:end-num_stim));
zauc = (crfresults.auc(1:end-num_stim,:) -mean(crfresults.auc(1:end-num_stim,:),1))./std(crfresults.auc(1:end-num_stim,:),1);

crfsigBK.zepsum = cat(1, crfsigBK.zepsum, zepsum);
crfsigBK.zauc = cat(1, crfsigBK.zauc, zauc);

if length(crfresults.epsum)-num_stim ~= Nneurons
    error('mismatch crfresults and expected number of neurons')
end

crfsigBK.ICencoder = cat(1, crfsigBK.ICencoder, ICencoder(neusigBK) );
crfsigBK.RCencoder = cat(1, crfsigBK.RCencoder, RCencoder(neusigBK) );
crfsigBK.inducerencoder = cat(1, crfsigBK.inducerencoder, inducerencoder(neusigBK) );

for g = 1:numel(neurongroupnames)
crfsigBK.(neurongroupnames{g}) = cat(1, crfsigBK.(neurongroupnames{g}), neurongroups.(neurongroupnames{g})(neusigBK) );
end

for isp =1:4
    switch isp
        case 1
            neuoi = true(Nneurons, 1);
            neutit = 'All';
        case 2
            neuoi = ICencoder(neusigBK);
            neutit = 'IC-enc';
        case 3
            neuoi = RCencoder(neusigBK);
            neutit = 'RC-enc';
        case 4
            neuoi = inducerencoder(neusigBK);
            neutit = 'ind-enc';
    end

Nneucum(ises,isp) = nnz(neuoi);
epsumacc(ises,isp) = mean(epsum(neuoi));
aucacc(ises,isp,:) = mean(auc(neuoi,:),1);
zepsumacc(ises,isp) = mean(zepsum(neuoi));
zaucacc(ises,isp,:) = mean(zauc(neuoi,:),1);
epmeanacc(ises,isp) = mean(epmean(neuoi));
end

end

%%
tempmat = epsumacc;

figure; plot(tempmat', 'o-')

temp = tempmat(all(~isnan(tempmat),2), :);
[P,ANOVATAB,STATS] = kruskalwallis(temp);
[c,m,h] = multcompare(STATS);
set(gca, 'YTick', h.CurrentAxes.YTick, 'YTickLabel', flip(ytl), 'FontSize', 11)
xlabel('Node Strength');
% set(h, 'Position',[300 200 400 100])

signrank(tempmat(:,1), tempmat(:,2))
signrank(tempmat(:,2), tempmat(:,3))
signrank(tempmat(:,2), tempmat(:,4))

%%
tempvec = crfsigBK.epmean;

figure; hold all
h = histogram(tempvec);
histogram(tempvec(crfsigBK.inducerencoder==1), 'BinEdges', h.BinEdges)
histogram(tempvec(crfsigBK.RCencoder==1), 'BinEdges', h.BinEdges)
histogram(tempvec(crfsigBK.ICencoder==1), 'BinEdges', h.BinEdges)

%%
pltopt = 3;
switch pltopt
    case 1
tempepagg = crfsigBK.epsum;
eplab = 'CRF Node Strength';
    case 2
tempepagg = crfsigBK.zepsum;
eplab = 'z-scored CRF Node Strength';
    case 3
tempepagg = crfsigBK.epmean;
eplab = 'CRF Mean Edge Potential';
end

Nneurons = length(crfsigBK.inducerencoder);
encgcols = [.5 .5 .5; 0 0.7 0; 1 0.5 0; 0 0 0];
Nneugcum = [];
neugcum = [];
epgcum = [];
aucgcum = [];
ytl = {};
for isp =1:4
    switch isp
        case 1
            neuoi = true(Nneurons, 1);
            neutit = 'sigBK';
        case 2
            neuoi = crfsigBK.ICencoder;
            neutit = 'IC-enc';
        case 3
            neuoi = crfsigBK.RCencoder;
            neutit = 'RC-enc';
        case 4
            neuoi = crfsigBK.inducerencoder;
            neutit = 'ind-enc';
    end
    neuoi = ismember(neusigBK, find(neuoi));
    if nnz(neuoi)==0
        continue
    end
    ytl = cat(1,ytl,neutit);
    Nneugcum = cat(1, Nneugcum, nnz(neuoi));
    neugcum = cat(1, neugcum, isp*ones(nnz(neuoi),1) );
    epgcum = cat(1, epgcum, tempepagg(neuoi));
    aucgcum = cat(1, aucgcum, crfsigBK.auc(neuoi,:));
end

[P,ANOVATAB,STATS] = kruskalwallis(epgcum, neugcum);
[c,m,h] = multcompare(STATS);
set(gca, 'YTick', h.CurrentAxes.YTick, 'YTickLabel', flip(ytl), 'FontSize', 11)
xlabel('Node Strength');
set(h, 'Position',[300 200 400 100])

fs=16;
figure('Position',[300*(isp-1), 100, 300 330])
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', sprintf('Aggregate N=%d Sessions', numel(ses2agg)), 'interpreter', 'none', 'edgecolor', 'none')
t = tiledlayout(1,1);
ax1 = axes(t);
xl = [0.5 length(ytl)+0.5];
hold on
b = boxplot(epgcum, neugcum, 'Colors', encgcols, 'notch', 'on', 'Symbol', 'k+');
set(b,{'linew'},{3})
yl = ylim;
text(xl(2), yl(1)+0.01*range(yl), sprintf('Pkw=%.4f', P), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fs)
xlim(xl)
% if doyl; ylim(yl); end
set(ax1, 'FontSize', fs, 'XTick', 1:length(ytl),  'XTickLabel', ytl, 'XTickLabelRotation', 45)
% xlabel(ax1, 'Neuron Subset')
ylabel(ax1, eplab)

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
