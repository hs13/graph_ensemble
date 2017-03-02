function [] = figS1_pa_whole_model_TF(param)

% parameters
expt_name = param.expt_name;
test_ee = param.test;
tf_seq = param.tf_seq;
ee = param.ee;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.pa;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
num_expt = length(expt_name);
linew = param.linew;
qnoise = param.qnoise;

load(ccode_path);
load(param.four_stim_cmap);

%% initialize
num_tf = length(tf_seq);
pred_mLL = cell(4,num_tf,num_expt);
pred_stats = zeros(2,num_expt,num_tf,3);
thr = zeros(num_expt,num_tf,1);
LL_pred = cell(4,num_tf,num_expt);

% cos_sim = {};
% cos_sim_avg = [];
% cos_thresh = [];
% pred = {};
% pred_stats = [];

%%
for n = 1:num_expt
    
    expt_ee = ee{n}{1};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\coords.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    
    % plot model
    if n==1
        coords = Coord_active;
        coords(end+1,:) = [0 max(coords(:,2))];
        coords(end+1,:) = [0 0];
        coords(end+1,:) = [max(coords(:,1)) 0];
        coords(end+1,:) = [max(coords(:,1)) max(coords(:,2))];
        figure; set(gcf,'color','w','position',[2084 472 402 319])
        plotGraphModel(best_model.graph,coords,best_model.edge_pot,[],gray(64))
        print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' ...
            expt_ee '_' ge_type '_model.pdf'])
    end
    
    for m = 1:num_tf
        
        load([data_path expt_name{n} '\' test_ee{n}{m} '.mat']);
        num_stim = length(setdiff(unique(vis_stim),0));
        data = data';
        vis_stim = vis_stim';
        num_frame = size(data,2);

        % calculate likelihood
        LL_frame = zeros(num_frame,num_stim);
        for ii = 1:num_frame
            for jj = 1:num_stim
                stim_vec = zeros(num_stim,1);
                stim_vec(jj) = 1;
                data_stim = [data(:,ii);stim_vec]';
                LL_frame(ii,jj) = compute_avg_log_likelihood(best_model.node_pot,...
                    best_model.edge_pot,best_model.logZ,data_stim);
            end
        end
        
        % make prediction
        [~,pred] = max(LL_frame,[],2);
        pred_final = zeros(size(pred));
        pred_mat = zeros(num_stim,num_frame);
        for ii = 1:num_stim % go through each stimulus
           pred_cr = pred==ii;
           for jj = setdiff(1:num_stim,ii) % 
                LLs = LL_frame(:,ii)-LL_frame(:,jj);
                % threshold by 3 std of noise
                th1 = quantile(LLs(:),qnoise);
                th2 = quantile(LLs(:),1-qnoise);
                LLs_th = LLs;
                LLs_th(LLs<=th1 | LLs>=th2) = NaN;
                thr = 3*nanstd(LLs_th(:))+nanmean(LLs_th(:));
                if ~isnan(thr)
                    pred_cr(LLs<thr) = 0;
                end
%                 LL_pred{ii,m,n} = LLs(vis_stim==ii);
            end
            pred_mat(ii,:) = pred_cr;
            pred_final = pred_final+pred_cr*ii;
            % collect LL per stim
            pred_mLL{ii,m,n} = LLs(vis_stim==ii);

        end
        
%         % ROC
%         auc = zeros(num_expt,2);
%         xx = cell(num_expt,2);
%         yy = cell(num_expt,2);
%         for ii = 1:num_expt
%             [auc(ii,1),xx{ii,1},yy{ii,1}] = plotROCmultic(LLs{ii,1},scores{ii,1},1,mycc.red_light,linew);
%         end
%         % calculate mean curve
%         xvec = 0:0.02:1;
%         ymat = zeros(num_expt,length(xvec),2);
%         for ii = 1:num_expt
%             [~,uid] = unique(xx{ii,1});
%             ymat(ii,:,1) = interp1(xx{ii,1}(uid),yy{ii,1}(uid),xvec);
%             [~,uid] = unique(xx{ii,2});
%             ymat(ii,:,2) = interp1(xx{ii,2}(uid),yy{ii,2}(uid),xvec);
%         end
        
        % prediction statistics
        for ii = 1:num_stim
            true_label = vis_stim'==ii;
            TP = sum(pred_final==ii & true_label==1);
            TN = sum(pred_final~=ii & true_label==0);
            FP = sum(pred_final==ii & true_label==0);
            FN = sum(pred_final~=ii & true_label==1);
            acc = (TP+TN)/(TP+TN+FN+FP);
            prc = TP/(TP+FP);
            rec = TP/(TP+FN);
            pred_stats(ii,n,m,:) = [acc,prc,rec];
        end

        % plot prediction
        if n==1
            plot_pred_raster(pred_mat,vis_stim',cmap);
            title(['TF = ' num2str(tf_seq(m))]);
            print(gcf,'-dpdf','-painters','-bestfit',[fig_path expt_name{n} '_' ...
                expt_ee '_' ge_type '_pred_TF' num2str(tf_seq(m)) '.pdf'])
        end
        
    end
    
end

%% plot stats
figure;
set(gcf,'color','w','position',[2041 533 993 235]);
set(gcf,'paperpositionmode','auto')

ww = 0.2;
stepsz = 0.2;

% accuracy
subplot(1,3,1); hold on
for m = 1:num_tf
    h = boxplot(reshape(pred_stats(:,:,m,1),[],1),'positions',m-stepsz,...
        'width',ww,'colors','k');
    setBoxStyle(h,linew)
end
xlim([0.5 m+0.5]); ylim([0 1])
set(gca,'xtick',1:m,'xticklabel',tf_seq);
xlabel('Temporal frequency (Hz)'); ylabel('Accuracy')
box off
% significance test
pval = zeros(1,num_tf);
for m = 2:num_tf
    pval(m) = ranksum(reshape(pred_stats(:,:,1,1),[],1),reshape(pred_stats(:,:,m,1),[],1));
end
title(num2str(pval(2:end)));

% precision
subplot(1,3,2); hold on
for m = 1:num_tf
    h = boxplot(reshape(pred_stats(:,:,m,2),[],1),'positions',m-stepsz,...
        'width',ww,'colors','k');
    setBoxStyle(h,linew)
end
xlim([0.5 m+0.5]); ylim([0 1])
set(gca,'xtick',1:m,'xticklabel',tf_seq);
xlabel('Temporal frequency (Hz)'); ylabel('Precision')
box off
% significance test
pval = zeros(1,num_tf);
for m = 2:num_tf
    pval(m) = ranksum(reshape(pred_stats(:,:,1,2),[],1),reshape(pred_stats(:,:,m,2),[],1));
end
title(num2str(pval(2:end)));

% recall
subplot(1,3,3); hold on
for m = 1:num_tf
    h = boxplot(reshape(pred_stats(:,:,m,3),[],1),'positions',m-stepsz,...
        'width',ww,'colors','k');
    setBoxStyle(h,linew)
end
xlim([0.5 m+0.5]); ylim([0 1])
set(gca,'xtick',1:m,'xticklabel',tf_seq);
xlabel('Temporal frequency (Hz)'); ylabel('Recall')
box off
% significance test
pval = zeros(1,num_tf);
for m = 2:num_tf
    pval(m) = ranksum(reshape(pred_stats(:,:,1,3),[],1),reshape(pred_stats(:,:,m,3),[],1));
end
title(num2str(pval(2:end)));

print(gcf,'-dpdf','-painters','-bestfit',[fig_path ge_type ...
    '_tf_pred_stats.pdf'])


end

