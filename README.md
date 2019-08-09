# DIRECT
Clustering on baseline phenotypes

1.[Pick number of archetypes](1.pick_archetypes_number.R) - pick few random seeds, and run the robustArchetypes method 10 times with each seed. Use screeplots and adjusted Rand index to pick optimal number of archetypes to use in clustering. 
  
  a. [Test archetype stability by archetype membership](1a.archetype_stability_by_threshold.R) - test clustering stability with adjusted Rand index at different thresholds of archetype membership

2.[Define 4 archetypes](2.archetypes_WP2.2.R) - plot PCA plot with colors representing the dominant archetype, shaded by archetype membership

3.[Find phenotypes correlated with archetypes](3.phenotype_differences.R) - heatmaps with differences

4.[Test progression differences](4.progression.R) - do archetypes correlate with HbA1c progression?

5.[Differences in genetics](5.psGRS.R) - test 6 partition T2D-GRSs

6.[Omics signature of archetypes](6.omics_signatures.R)

  a. [Select omics variables to plot](6a.select_top_omics_top_plot.sh)

  b. [Heatmaps of omics signatures](6b.omics_heatmap.R)

7.[Immune cell signatures in DEGs](7.immune_cell_enrichment.R)

8.[Longitudinal archetypes analysis](8.archetypes_longitudinal.R)
