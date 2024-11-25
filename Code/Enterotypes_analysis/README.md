# Enterotype Analysis Code

In the `Code/ISN_construction/` folder, there are the early stages of the Enterotype analysis, i.e., the calculation of the population-based networks and the ISNs that are the input for the analysis, referred to in the "Enterotype-based analysis" in the paper.  
The analysis in this folder, on the contrary, are post-processing analyses after the Encoding/Decoding pipeline (PLEX.I) has been run. We refer to the Plex.I paper (https://pubmed.ncbi.nlm.nih.gov/37928248/) for details on the Encoding/Decoding procedures. 

## Input
- Mapping of node name to the OTU.
- Taxonomic annotation: taxa → genus → family → … 
- Encoding/Decoding results: A dataframe with three columns: node name, the rank of the sum of the distances involving the node, and the p-value.

## Aim
- Identify key taxa that have significantly different connections (i.e., edges) before/after the treatment, divided based on the trajectory of individuals' Enterotypes.

## Structure

- `Enterotype_Analysis.r` details all the steps.
- To help the reader understand the various steps, we generate a Markdown file showcasing the intermediate and final results: `Enterotype_Analysis.html` and `Enterotype_Analysis.rmd`.
