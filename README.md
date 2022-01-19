# SpatialClustersBreastCancerIndividualData


This repository contains the code used for the analysis in "A Flexible Method for Identifying Spatial Clusters of Breast Cancer Using Individual-Level Data" by Kamenetsky, Trentham-Dietz, Newcomb, Zhu, Gangnon (2022).

The files used in the analysis are described below:

- *0_descriptivetables.R*: Creates Table 1 (along with other tables).
- *0_exploratoryviz.R*: Exploratory data analysis. Creates Figure 2.
- *1_imputation.R*: Imputation of missing data using predictive mean
  matching.
- *1_nbdanalysis.R*: Spatial cluster analysis by each Wisconsin
  county-neighborhood.
- *2_commutingcalcs.R*: Calculates median distance traveled in Wisconsin using
  data from the National Household Travey Survey, 2017.
- *3_logisticmodelsimputed_bic.R*: Performs the individual-level modeling with
  identified clusters across fully-adjusted and age-adjusted model.
- *4_clustermap.R*: Several visualizations including Figures 1 and 3.
- *4_results_clustersidentify.R*: Identifies the clusters from the model and
  creates `clusterdf_wisc_bic.csv` which is used in subsequent figures and
analyses.
- *4_resultstable_uniqwomen.R*: Creates Table 2.
- *5_resultsviz_micromaps_serverplotgrid.R*: Creates Figure 4.




The data comes from the Wisconsin Women's Health Study. Due to privacy issues,
it is not publicly available.

Code developed by and repository maintained by M. Kamenetsky.
