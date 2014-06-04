Analysis for Fig 3, S5 and Supplemental Table 1

To carry out k-means and associated analyses, the scripts in this folder were run in the following order and serve the following functions (even more details can be found at the top of the script files):

1. k-means_data_prep.m:
	-averages the replicates and stores to new data structures
	-has functionality to plot abundance by location for all 
	proteins, or proteins between given thresholds

2. k_means.m:
	-runs k-means clustering on protein abundance across locations
	-filters out proteins in identified low-abundance clusters
	-re-runs k-means on high abundance proteins normalized by zscore

3. finding_sqeucdistances.m:
	-finds square Euclidean distance between clusters of two given 
	colonization states
	-reports protein IDs for entered cluster number, and re-plots 
	abundance over location from original data as a sanity check
	-has functionality to identify protein overlap for two clusters 
	entered from two different colonization states
	


NOTE: Because the cluster numbers change each time the script is run, the data generated for the final figures from k-means clustering is saved in the following .mat file:

clustering_work_with_final_norm.mat