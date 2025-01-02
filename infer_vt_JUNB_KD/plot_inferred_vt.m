%% visualize inferred transcription rates for three classes of genes: DDGs, delayed and biphasic
% Lists of gene names in 3 classes:
% DDGs: Delayed_genes_high_FC6.csv
% Delayed: rejected_delayed_D2D.csv
% Biphasic: rejected_biphasic_D2D.csv

close all
clear

% load inference results 
result_name = {'results_ddgs.mat','results_delay.mat','results_biphasic.mat'};

% 3 classes of genes
ftitle = {'DDGs','Delayed','Biphasic'};

% calculate the confidence intervals of inferred transcription rates and
% visuaize the final results. Parallel computing toolbox was used.
for ii = 1:3

    % load inference results
    load(result_name{ii})
    
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
    
    % time vector
    t = 0:15:720;

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