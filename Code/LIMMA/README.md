# LIMMA Calculation Code Documentation

## Overview
This script is designed to perform differential network analysis using the LIMMA (Linear Models for Microarray Data) method, and is referred in section "Differential network analysis" of the paper. It calculates and analyzes the differences in microbial network connectivity between two conditions (responders and non-responders) based on specific outcomes (e.g., clinical, biomarker, or endoscopic outcomes). It processes the microbiome data, identifies significant edges (connections between taxa), and visualizes the results.

## Input
- **Taxonomic annotation**: Metadata detailing taxa → genus → family → (higher taxonomic levels).
- **ISN network**: Key input, for each iteration, is the ISN (outcome of the code in the ISN_construction folder).
- **Metadata file (`w=_n355_metadata.txt`)**: Metadata files that contain information about each individual sample, including their outcomes and response to treatment. These are used to group the samples into responders and non-responders.

## Goal
The primary goal of this analysis is to compute and identify significant microbial network edges that differ between two groups (responders vs. non-responders). The LIMMA method is applied to calculate differential expression (logFC) for each edge in the network, followed by statistical significance testing (adjusted for multiple testing). The results are visualized and saved for further interpretation.

## Structure of the Code

- Limma_calculation.r Key Components

> Import Libraries, file and metadata loading; yhe necessary R libraries are loaded to perform the analysis:

lionessR: For network analysis.
igraph: For network graph construction and visualization.
reshape2, dplyr: For data manipulation.
limma: For performing linear modeling and statistical analysis.


> Differential Network Calculation (LIMMA)

The script calculates the difference in network connectivity between two groups (responders vs. non-responders) for each edge. The edges with an absolute difference greater than a threshold (0.25, but we also tested 0.5) are considered significant.
The LIMMA method is then applied to the filtered network to identify which edges are differentially expressed between the two groups, adjusting for multiple testing using FDR (False Discovery Rate) correction. Steps are specular to https://github.com/mararie/lionessR 

> Visualization
The significant network edges are visualized using igraph. The visualization highlights the most significant edges and generates a graphical representation of the network.
Networks are plotted for different sets of edges:
All edges
Top 50 edges based on logFC
Significant edges with p-value < 0.05

## Output

The results are saved, for each iteartion, in several formats in Results/LIMMA and Graphs/LIMMA:
A file with all edges with the corresponding p-values.
A file with only the top 50 edges.
A file with only the significant edges (p-value < 0.05).
The network plots are saved as PNG files for visual inspection.

Finally, Limma_creating_LIMMA_folders_moving_there_all_step2.r moves all the results in a single, common, folder to help the postprocessing (i.e. creating a single file comparing the LIMMA results between different cohorts and treatments). 

An example of the result is shown here, after the postprocessing and aggregation done with the scripts in Code/Grouping_LIMMA_results_shared_taxon for LIMMA with 0.25 as the threshold, at week 0 and for the endoscopic outcome. The taxa are aggregated at the family level.
![image](../../images/Limma_025_family_endoscopic_w0.png)

