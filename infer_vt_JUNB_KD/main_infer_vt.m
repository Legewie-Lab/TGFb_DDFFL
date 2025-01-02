%% main script for transcription rate inference
close all
clear

% import RNA-seq data of JUNB KD; t1:1.5h; t3:3h; t4:t6h; t2:12h
data_kd_t1 = readtable('KD_H_t1vsKD_H_ctrl_deg.csv');
data_kd_t2 = readtable('KD_H_t2vsKD_H_ctrl_deg.csv');
data_kd_t3 = readtable('KD_H_t3vsKD_H_ctrl_deg.csv');
data_kd_t4 = readtable('KD_H_t4vsKD_H_ctrl_deg.csv');

% import RNA-seq data of JUNB KD control; t1:1.5h; t3:3h; t4:t6h; t2:12h
data_ct_t1 = readtable('KD_N_t1vsKD_N_ctrl_deg.csv');
data_ct_t2 = readtable('KD_N_t2vsKD_N_ctrl_deg.csv');
data_ct_t3 = readtable('KD_N_t3vsKD_N_ctrl_deg.csv');
data_ct_t4 = readtable('KD_N_t4vsKD_N_ctrl_deg.csv');

% get common genes recorded across all time points in JUNB KD and KD control conditions 
gene_id_unique = unique([data_kd_t1.gene_id;data_kd_t2.gene_id;data_kd_t3.gene_id;data_kd_t4.gene_id;...
                         data_ct_t1.gene_id;data_ct_t2.gene_id;data_ct_t3.gene_id;data_ct_t4.gene_id]);

ind_kd_t1 = ismember(gene_id_unique,data_kd_t1.gene_id);
ind_kd_t2 = ismember(gene_id_unique,data_kd_t2.gene_id);
ind_kd_t3 = ismember(gene_id_unique,data_kd_t3.gene_id);
ind_kd_t4 = ismember(gene_id_unique,data_kd_t4.gene_id);

ind_ct_t1 = ismember(gene_id_unique,data_ct_t1.gene_id);
ind_ct_t2 = ismember(gene_id_unique,data_ct_t2.gene_id);
ind_ct_t3 = ismember(gene_id_unique,data_ct_t3.gene_id);
ind_ct_t4 = ismember(gene_id_unique,data_ct_t4.gene_id);

ind_common = ind_kd_t1 & ind_kd_t2 & ind_kd_t3 & ind_kd_t4 ...
           & ind_ct_t1 & ind_ct_t2 & ind_ct_t3 & ind_ct_t4;

gene_id_common = gene_id_unique(ind_common);

% filtering and sorting data by the common genes
data_kd_t1 = sortrows(data_kd_t1(ismember(data_kd_t1.gene_id,gene_id_common),:),1);
data_kd_t2 = sortrows(data_kd_t2(ismember(data_kd_t2.gene_id,gene_id_common),:),1);
data_kd_t3 = sortrows(data_kd_t3(ismember(data_kd_t3.gene_id,gene_id_common),:),1);
data_kd_t4 = sortrows(data_kd_t4(ismember(data_kd_t4.gene_id,gene_id_common),:),1);

data_ct_t1 = sortrows(data_ct_t1(ismember(data_ct_t1.gene_id,gene_id_common),:),1);
data_ct_t2 = sortrows(data_ct_t2(ismember(data_ct_t2.gene_id,gene_id_common),:),1);
data_ct_t3 = sortrows(data_ct_t3(ismember(data_ct_t3.gene_id,gene_id_common),:),1);
data_ct_t4 = sortrows(data_ct_t4(ismember(data_ct_t4.gene_id,gene_id_common),:),1);

gene_id = data_kd_t1.gene_id;

% measured time points
tseq = [1.5 3 6 12]*60;

% obtain log2 fold change in KD and control conditions
gene_exp_kd = [data_kd_t1.log2FoldChange,...
               data_kd_t3.log2FoldChange,...
               data_kd_t4.log2FoldChange,...
               data_kd_t2.log2FoldChange];

gene_exp_ct = [data_ct_t1.log2FoldChange,...
               data_ct_t3.log2FoldChange,...
               data_ct_t4.log2FoldChange,...
               data_ct_t2.log2FoldChange];

n_genes = length(gene_id); % number of common genes

% time points for inference
t = 0:15:720;

% create rna_seq structure to store all data needed
rna_seq = struct;
rna_seq.t = tseq;

rna_seq.fold_change_ct = 2.^gene_exp_ct;
rna_seq.fold_change_kd = 2.^gene_exp_kd;

% pre-allocate variables for x (mRNA fold change) and vt (transcription rate)
x_ct = zeros(n_genes,length(t));
x_kd = zeros(n_genes,length(t));

vt_ct = zeros(n_genes,length(t));
vt_kd = zeros(n_genes,length(t));

% use the mean mRNA lifetime 9h 
tau = 9*60;

% inference for each gene
parfor ii = 1:n_genes
    
    [x_ct(ii,:),vt_ct(ii,:),x_kd(ii,:),vt_kd(ii,:)] = calc_x_t(ii,tau,rna_seq,t);

end

% store the inference results in the rna_seq structure
rna_seq.x_ct = x_ct;
rna_seq.vt_ct = vt_ct;

rna_seq.x_kd = x_kd;
rna_seq.vt_kd = vt_kd;

% save inference results
result_name = 'infer_vt.mat';
save(result_name)

%% visualize inferred transcription rates for three classes of genes: DDGs, delayed and biphasic
% Lists of gene names in 3 classes:
% DDGs: Delayed_genes_high_FC6.csv
% Delayed: rejected_delayed_D2D.csv
% Biphasic: rejected_biphasic_D2D.csv

% 3 classes of genes
ftitle = {'DDGs','Delayed','Biphasic'};

% calculate the confidence intervals of inferred transcription rates and
% visuaize the final results. Parallel computing toolbox was used.
for ii = 1:3
   
    if ii == 1
        load(result_name)
        % names of all inferred genes
        gene_name = data_ct_t1.gene_name;
        % names of DDGs which is listed in the file 'Delayed_genes_high_FC6.csv'
        data_ddg = readtable('Delayed_genes_high_FC6.csv');
        % names DDGs overlapping with all inferred genes 
        [gene_name_ddg,ind_ddg,~] = intersect(gene_name,data_ddg.gene_name);
        % get the inference results and data for DDGs 
        rna_seq.fold_change_ct=rna_seq.fold_change_ct(ind_ddg,:);
        rna_seq.fold_change_kd=rna_seq.fold_change_kd(ind_ddg,:);
        rna_seq.x_ct=rna_seq.x_ct(ind_ddg,:);
        rna_seq.x_kd=rna_seq.x_kd(ind_ddg,:);
        rna_seq.vt_ct=rna_seq.vt_ct(ind_ddg,:);
        rna_seq.vt_kd=rna_seq.vt_kd(ind_ddg,:);
        gene_id = gene_id(ind_ddg);
        gene_name = gene_name(ind_ddg);

    elseif ii == 2
        load(result_name)
        % names of all inferred genes
        gene_name = data_ct_t1.gene_name;
        % names of delayed genes which is listed in the file 'rejected_delayed_D2D.csv'
        data_delay = readtable('rejected_delayed_D2D.csv');
        % names delayed genes overlapping with all inferred genes
        [gene_name_delay,ind_delay,~] = intersect(gene_name,data_delay.gene_name);
        % get the inference results and data for delayed genes
        rna_seq.fold_change_ct=rna_seq.fold_change_ct(ind_delay,:);
        rna_seq.fold_change_kd=rna_seq.fold_change_kd(ind_delay,:);
        rna_seq.x_ct=rna_seq.x_ct(ind_delay,:);
        rna_seq.x_kd=rna_seq.x_kd(ind_delay,:);
        rna_seq.vt_ct=rna_seq.vt_ct(ind_delay,:);
        rna_seq.vt_kd=rna_seq.vt_kd(ind_delay,:);
        gene_id = gene_id(ind_delay);
        gene_name = gene_name(ind_delay);

    elseif ii == 3
        load(result_name)
        % names of all inferred genes
        gene_name = data_ct_t1.gene_name;
        % names of biphasic genes which is listed in the file 'rejected_biphasic_D2D.csv'
        data_biph = readtable('rejected_biphasic_D2D.csv');	
        % names biphasic genes overlapping with all inferred genes
        [gene_name_biph,ind_biph,~] = intersect(gene_name,data_biph.gene_name);
        % get the inference results and data for biphasic genes
        rna_seq.fold_change_ct=rna_seq.fold_change_ct(ind_biph,:);
        rna_seq.fold_change_kd=rna_seq.fold_change_kd(ind_biph,:);
        rna_seq.x_ct=rna_seq.x_ct(ind_biph,:);
        rna_seq.x_kd=rna_seq.x_kd(ind_biph,:);
        rna_seq.vt_ct=rna_seq.vt_ct(ind_biph,:);
        rna_seq.vt_kd=rna_seq.vt_kd(ind_biph,:);
        gene_id = gene_id(ind_biph);
        gene_name = gene_name(ind_biph);
    end
    
    % separate up- and down-regulated genes
    ind_up_reg = rna_seq.fold_change_ct(:,end) > 1;
    ind_dn_reg = rna_seq.fold_change_ct(:,end) <= 1;
    
    vt_ct_up = rna_seq.vt_ct(ind_up_reg,:);
    vt_ct_dn = rna_seq.vt_ct(ind_dn_reg,:);

    vt_kd_up = rna_seq.vt_kd(ind_up_reg,:);
    vt_kd_dn = rna_seq.vt_kd(ind_dn_reg,:);

    % calculate the confidence intervals for inferred transcription rates
    % via bootstrapping (2000 boostrap samples)
    nboot = 2000;
    options=statset(UseParallel=true);

    cal_fun = @(x) median(x);

    [ci_vt_ct_up{ii},bootstat_vt_ct_up{ii}] = bootci(nboot,{cal_fun,vt_ct_up},'Options',options);
    [ci_vt_ct_dn{ii},bootstat_vt_ct_dn{ii}] = bootci(nboot,{cal_fun,vt_ct_dn},'Options',options);
    [ci_vt_kd_up{ii},bootstat_vt_kd_up{ii}] = bootci(nboot,{cal_fun,vt_kd_up},'Options',options);
    [ci_vt_kd_dn{ii},bootstat_vt_kd_dn{ii}] = bootci(nboot,{cal_fun,vt_kd_dn},'Options',options);

    
    [ci_vt_diff_up{ii},bootstat_vt_diff_up{ii}] = bootci(nboot,{cal_fun,vt_kd_up-vt_ct_up},'Options',options);
    [ci_vt_diff_dn{ii},bootstat_vt_diff_dn{ii}] = bootci(nboot,{cal_fun,vt_kd_dn-vt_ct_dn},'Options',options);

    % visualizing results
    fa = 0.3; % alpha controling transparency of the confidence bands
    gc = [0 200 0]/255; % color
    
    % up-regulated genes
    figure(1)
    subplot(3,1,ii)
        yline(0,'Color',[1 1 1]*0.75,'LineStyle','--')
    hold on
    fill([t,fliplr(t)],[ci_vt_diff_up{ii}(1,:),fliplr(ci_vt_diff_up{ii}(2,:))],...
         gc,'LineStyle','none','FaceAlpha',fa)
    plot(t,cal_fun(vt_kd_up-vt_ct_up),'color',gc)
    xlim([-0.02 720])
    xticks(0:120:720)
    ylim([-4.02 2.02])
    yticks(-4:2:2)
    set(gca,'tickdir','out','box','off')
    xlabel('Time (min)')
    title(ftitle{ii},'FontWeight','normal')

    % down-reglated genes
    figure(2)
    subplot(3,1,ii);
    yline(0,'Color',[1 1 1]*0.75,'LineStyle','--')
    hold on
    fill([t,fliplr(t)],[ci_vt_diff_dn{ii}(1,:),fliplr(ci_vt_diff_dn{ii}(2,:))],...
         gc,'LineStyle','none','FaceAlpha',fa)
    plot(t,cal_fun(vt_kd_dn-vt_ct_dn),'color',gc)
    xlim([-0.02 720])
    xticks(0:120:720)
    ylim([-0.5 1])
    yticks(-0.5:0.5:1)
    set(gca,'tickdir','out','box','off')
    xlabel('Time (min)')
    title(ftitle{ii},'FontWeight','normal')

end

% setting figure size
figure(1)
width=4;
height=10.5;
fig_position = [0,0,width,height];
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'position',fig_position,'color','white');
set(gcf,'paperposition',fig_position,'PaperSize',fig_position(3:4));

figure(2)
width=4;
height=10.5;
fig_position = [0,0,width,height];
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'position',fig_position,'color','white');
set(gcf,'paperposition',fig_position,'PaperSize',fig_position(3:4));


%%  
% The function that performs the inference
function [x_ct,vt_ct,x_kd,vt_kd] = calc_x_t(ii,tau,rna_seq,t)

    % calculate piecewise splines using Modified Akima Interpolation for control condition 
    spline_ct = makima([0 rna_seq.t],[1 rna_seq.fold_change_ct(ii,:)]);
    x_ct = max([zeros(size(ppval(spline_ct,t)));ppval(spline_ct,t)]);

    % calculate piecewise splines using Modified Akima Interpolation for KD condition 
    spline_kd = makima([0 rna_seq.t],[1 rna_seq.fold_change_kd(ii,:)]);
    x_kd = max([zeros(size(ppval(spline_kd,t)));ppval(spline_kd,t)]);  

    % inference of transcription rate of control and KD 
    vt_ct = [1, diff(x_ct)./diff(t) * tau + x_ct(1:end-1)];
    vt_kd = [1, diff(x_kd)./diff(t) * tau + x_kd(1:end-1)];
    
    % avoid negative rates
    vt_ct(vt_ct<=0) = 1e-6;
    vt_kd(vt_kd<=0) = 1e-6;
  
end