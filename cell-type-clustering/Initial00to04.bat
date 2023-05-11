:: Automates running the first five scripts. Make sure all scripts are updated with the correct dataset name, and dataset.R is also updated. Make sure 03 and 04 are set at the whole level.

Rscript 00_convert_to_rds.R
Rscript 01_plot_qc.R
Rscript 02_compute_marker_scores.R
Rscript 03_analyze_barseq.R
Rscript 04_plot_analysis.R