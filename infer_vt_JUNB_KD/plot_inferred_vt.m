%%
close all
clear

result_name = 'infer_vt.mat';

ftitle = {'DDGs','Delayed','Biphasic'};

for ii = 1:3

   
    if ii == 1
        load(result_name)
        gene_name = data_ct_t1.gene_name;
        data_ddg = readtable('Delayed_genes_high_FC6.csv');
        [gene_name_ddg,ind_ddg,~] = intersect(gene_name,data_ddg.gene_name);
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
        gene_name = data_ct_t1.gene_name;

        data_delay = readtable('rejected_delayed_D2D.csv');
        [gene_name_delay,ind_delay,~] = intersect(gene_name,data_delay.gene_name);
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
        gene_name = data_ct_t1.gene_name;
        data_biph = readtable('rejected_biphasic_D2D.csv');	
        [gene_name_biph,ind_biph,~] = intersect(gene_name,data_biph.gene_name);
        rna_seq.fold_change_ct=rna_seq.fold_change_ct(ind_biph,:);
        rna_seq.fold_change_kd=rna_seq.fold_change_kd(ind_biph,:);
        rna_seq.x_ct=rna_seq.x_ct(ind_biph,:);
        rna_seq.x_kd=rna_seq.x_kd(ind_biph,:);
        rna_seq.vt_ct=rna_seq.vt_ct(ind_biph,:);
        rna_seq.vt_kd=rna_seq.vt_kd(ind_biph,:);
        gene_id = gene_id(ind_biph);
        gene_name = gene_name(ind_biph);
    end
    
    ind_up_reg = rna_seq.fold_change_ct(:,end) > 1;
    ind_dn_reg = rna_seq.fold_change_ct(:,end) <= 1;

    vt_ct_up = rna_seq.vt_ct(ind_up_reg,:);
    vt_ct_dn = rna_seq.vt_ct(ind_dn_reg,:);

    vt_kd_up = rna_seq.vt_kd(ind_up_reg,:);
    vt_kd_dn = rna_seq.vt_kd(ind_dn_reg,:);

    nboot = 2000;
    options=statset(UseParallel=true);

    cal_fun = @(x) median(x);

    [ci_vt_ct_up{ii},bootstat_vt_ct_up{ii}] = bootci(nboot,{cal_fun,vt_ct_up},'Options',options);
    [ci_vt_ct_dn{ii},bootstat_vt_ct_dn{ii}] = bootci(nboot,{cal_fun,vt_ct_dn},'Options',options);
    [ci_vt_kd_up{ii},bootstat_vt_kd_up{ii}] = bootci(nboot,{cal_fun,vt_kd_up},'Options',options);
    [ci_vt_kd_dn{ii},bootstat_vt_kd_dn{ii}] = bootci(nboot,{cal_fun,vt_kd_dn},'Options',options);

    
    [ci_vt_diff_up{ii},bootstat_vt_diff_up{ii}] = bootci(nboot,{cal_fun,vt_kd_up-vt_ct_up},'Options',options);
    [ci_vt_diff_dn{ii},bootstat_vt_diff_dn{ii}] = bootci(nboot,{cal_fun,vt_kd_dn-vt_ct_dn},'Options',options);

    
    fa = 0.3;
    gc = [0 200 0]/255;

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