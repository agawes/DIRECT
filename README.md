# DIRECT
Clustering on baseline phenotypes

1.[Pick number of archetypes](1.pick_archetypes_number.R) - pick few random seeds, and run the robustArchetypes method 10 times with each seed. Use screeplots and adjusted Rand index to pick optimal number of archetypes to use in clustering. 

2.[Define 4 archetypes](2.archetypes_WP2.2.R) - plot PCA plot with colors representing the dominant archetype, shaded by archetype membership

3.[Find phenotypes correlated with archetypes](3.phenotype_differences.R) - heatmaps with differences

4.[Test progression differences](progression.R) - do archetypes correlate with HbA1c progression?

5.[Differences in genetics](psGRS.R) - test 6 partition T2D-GRSs

6.[Omics signature of archetypes]()

  a. [Select omics variables to plot]()

  b. [Heatmaps of omics signatures](omics_heatmap.R)

7.[Immune cell signatures in DEGs](immune_cell_enrichment.R)

8.[Longitudinal archetypes analysis](archetypes_longitudinal.R)
