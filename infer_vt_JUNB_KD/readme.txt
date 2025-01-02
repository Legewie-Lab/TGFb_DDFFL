This directory contains the Matlab scripts (main_infer_vt.m and plot_inferred_vt.m) and datasets used for inference of the input transcription rates of genes in three classes, i.e., DDGs, delayed and biphasic, upon JUNB knockdown (KD). The Parallel Computing Toolbox was used in the Matlab script to speedup the calculation. 

Details of data files:
Matlab scripts and data
main_infer_vt.m - main script for computational inference of input transcription rate.
plot_inferred_vt.m - calculate the confidence bands for the median of the inferred transcription rate difference (KD - control) in each gene class and visualize the final results.
results_ddgs.mat - inference result of time-dependent transcription rates (JUNB KD and control) for DDGs.
results_delay.mat - inference result of time-dependent transcription rates (JUNB KD and control) for delayed genes.
results_biphasic.mat - inference result of time-dependent transcription rates (JUNB KD and control) for biphasic genes.

RNA-seq data for JUNB KD control at different time points
KD_N_t1vsKD_N_ctrl_deg.csv - 1.5h
KD_N_t2vsKD_N_ctrl_deg.csv - 12h
KD_N_t3vsKD_N_ctrl_deg.csv - 3h
KD_N_t4vsKD_N_ctrl_deg.csv - 6h

RNA-seq data for JUNB KD at different time points
KD_H_t1vsKD_H_ctrl_deg.csv - 1.5h
KD_H_t2vsKD_H_ctrl_deg.csv - 12h
KD_H_t3vsKD_H_ctrl_deg.csv - 3h
KD_H_t4vsKD_H_ctrl_deg.csv - 6h

Gene lists
Delayed_genes_high_FC6.csv - list of DDGs
rejected_delayed_D2D.csv - list of delayed genes
rejected_biphasic_D2D.csv - list of biphasic genes
