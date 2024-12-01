%%
close all
clear

data_kd_t1 = readtable('KD_H_t1vsKD_H_ctrl_deg.csv');
data_kd_t2 = readtable('KD_H_t2vsKD_H_ctrl_deg.csv');
data_kd_t3 = readtable('KD_H_t3vsKD_H_ctrl_deg.csv');
data_kd_t4 = readtable('KD_H_t4vsKD_H_ctrl_deg.csv');

data_ct_t1 = readtable('KD_N_t1vsKD_N_ctrl_deg.csv');
data_ct_t2 = readtable('KD_N_t2vsKD_N_ctrl_deg.csv');
data_ct_t3 = readtable('KD_N_t3vsKD_N_ctrl_deg.csv');
data_ct_t4 = readtable('KD_N_t4vsKD_N_ctrl_deg.csv');

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

data_kd_t1 = sortrows(data_kd_t1(ismember(data_kd_t1.gene_id,gene_id_common),:),1);
data_kd_t2 = sortrows(data_kd_t2(ismember(data_kd_t2.gene_id,gene_id_common),:),1);
data_kd_t3 = sortrows(data_kd_t3(ismember(data_kd_t3.gene_id,gene_id_common),:),1);
data_kd_t4 = sortrows(data_kd_t4(ismember(data_kd_t4.gene_id,gene_id_common),:),1);

data_ct_t1 = sortrows(data_ct_t1(ismember(data_ct_t1.gene_id,gene_id_common),:),1);
data_ct_t2 = sortrows(data_ct_t2(ismember(data_ct_t2.gene_id,gene_id_common),:),1);
data_ct_t3 = sortrows(data_ct_t3(ismember(data_ct_t3.gene_id,gene_id_common),:),1);
data_ct_t4 = sortrows(data_ct_t4(ismember(data_ct_t4.gene_id,gene_id_common),:),1);

gene_id = data_kd_t1.gene_id;

tseq = [1.5 3 6 12]*60;

gene_exp_kd = [data_kd_t1.log2FoldChange,...
               data_kd_t3.log2FoldChange,...
               data_kd_t4.log2FoldChange,...
               data_kd_t2.log2FoldChange];

gene_exp_ct = [data_ct_t1.log2FoldChange,...
               data_ct_t3.log2FoldChange,...
               data_ct_t4.log2FoldChange,...
               data_ct_t2.log2FoldChange];

n_genes = length(gene_id);

t = 0:15:720;

rna_seq = struct;

rna_seq.t = tseq;

rna_seq.fold_change_ct = 2.^gene_exp_ct;
rna_seq.fold_change_kd = 2.^gene_exp_kd;

x_ct = zeros(n_genes,length(t));
x_kd = zeros(n_genes,length(t));

vt_ct = zeros(n_genes,length(t));
vt_kd = zeros(n_genes,length(t));

tau = 9*60;

parfor ii = 1:n_genes
    
    [x_ct(ii,:),vt_ct(ii,:),x_kd(ii,:),vt_kd(ii,:)] = calc_x_t(ii,tau,rna_seq,t);

end

rna_seq.x_ct = x_ct;
rna_seq.vt_ct = vt_ct;

rna_seq.x_kd = x_kd;
rna_seq.vt_kd = vt_kd;

save('infer_vt')


%%  
function [x_ct,vt_ct,x_kd,vt_kd] = calc_x_t(ii,tau,rna_seq,t)

    spline_ct = makima([0 rna_seq.t],[1 rna_seq.fold_change_ct(ii,:)]);
    x_ct = max([zeros(size(ppval(spline_ct,t)));ppval(spline_ct,t)]);  

    spline_kd = makima([0 rna_seq.t],[1 rna_seq.fold_change_kd(ii,:)]);
    x_kd = max([zeros(size(ppval(spline_kd,t)));ppval(spline_kd,t)]);  

    vt_ct = [1, diff(x_ct)./diff(t) * tau + x_ct(1:end-1)];
    vt_kd = [1, diff(x_kd)./diff(t) * tau + x_kd(1:end-1)];
    
    vt_ct(vt_ct<=0) = 1e-6;
    vt_kd(vt_kd<=0) = 1e-6;
  
end