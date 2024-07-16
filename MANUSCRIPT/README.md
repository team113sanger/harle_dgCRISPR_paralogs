# README

The `PLOTS` directory contatins plots used in the manuscript which are generated using R scripts in the `SCRIPTS` directory.

## Script descriptions and relationship to figures

This section describes the purpose of all scripts in the `SCRIPTS` directory and their respective manuscript figures.

The scripts source `prepare_data.R` to provide common paths, palettes, input data sets and helper functions. Input data are combined data sets which can be found in the `DATA/postprocessing` directory.

* **Figure 2** 
  * Script: `MANUSCRIPT/SCRIPTS/F2__heatmap_of_hits_with_gi_and_waterfall.R`
  * Output(s):
    * `MANUSCRIPT/PLOTS/F2__heatmap_without_total_hits_portrait.png`
  * Note: this script can also be used to generate a landscape version of the plots for presentations
  
* **Figure 4**
  * Script(s):
    * A: `MANUSCRIPT/SCRIPTS/F4A__gi_of_hits.R`
    * B: `MANUSCRIPT/SCRIPTS/F4B__pie_all_gene_pairs.R`
    * C: `MANUSCRIPT/SCRIPTS/F4C__number_of_times_a_hit_pie.R`
    * D: `MANUSCRIPT/SCRIPTS/F4D__pie_chart_of_hits_by_cancer_type.R`
  * Output(s):
    * A: `F4A__all_gi_scores.png`
    * B: `F4B__pie_chart_frequency_of_all_pairs.png`
    * C: `F4C__pie_chart_hit_frequency_altered_colours.png`
    * D: `MANUSCRIPT/PLOTS/F4D__pie_chart_of_hits_by_cancer_type.png`

* **Figure 5**
  * Script:
    * A: `MANUSCRIPT/SCRIPTS/F5A__correlation_of_median_hit_GI_with_number_of_hits.R`
    * B: `MANUSCRIPT/SCRIPTS/F5B__range_of_GI_scores.R`
    * C: `MANUSCRIPT/SCRIPTS/F5C__ranking_of_top_50percent_hits.R`
    * D: 
  * Output(s):
    * A: `F5A__median_GI_of_hits_per_pair_v_total_hits.png`
    * B: `F5B__GI_range_of_hits_per_cell_line.png`
    * C: `F5C__cell_line_ranking_of_hits_in_more_than_50pct_of_lines.png`
    * D:

* **Figure 6**
  * Script:
    * 
  * Output(s):
    * 
  * 

* **Figure 7**
  * Script: ``
  * Output(s):
    * 
  * 

* **Extended Data Figure 3**
  * Script(s): 
    * A: `SCRIPTS/QC/03_control_samples.R`
    * B: `SCRIPTS/QC/05_sample_essentiality.R`
    * C and D: `SCRIPTS/QC/02_library_statistics.R`
  * Output(s):
    * A: `MANUSCRIPT/PLOTS/total_read_pairs_per_sample.png`
    * B: `MANUSCRIPT/PLOTS/normalised_LFC_NNMD.png`
    * C: `MANUSCRIPT/PLOTS/median_normalised_counts_per_sample.png`
    * D: `MANUSCRIPT/PLOTS/gini_index_per_sample.png`

* **Extended Data Figure 4**
  * Script: `SCRIPTS/QC/03_control_samples.R`
  * Output(s):
    * `MANUSCRIPT/PLOTS/control_correlation_normalised_counts.png`

* **Extended Data Figure 5**
  * Script: 
    * A: `SCRIPTS/postprocessing/01_download_depmap_relative_copy_number.R`
    * B: `SCRIPTS/postprocessing/02_download_cmp_total_copy_number.R`
  * Output(s):
    * A: `MANUSCRIPT/PLOTS/DepMap_relative_copy_number.png`
    * B: `MANUSCRIPT/PLOTS/CMP_total_copy_number.png`

* **Extended Data Figure 6**
* Script: 
    * B: ``
  * Output(s):
    * B: ``
