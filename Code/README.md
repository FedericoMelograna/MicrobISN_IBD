# Code
 Code, for the ISN pipeline in Microbiome IBD data. In this folder, we collected the code to reproduce the various analysis. 

> In the Alpha_beta diversity subfolder, there is the code to analyse Alpha and Beta diversity for the vairous combination of cohorts + time + treatment. Moreover, in the same folder, there are also the exploratory analysis on the abundance/fractional abundance, both represented as a stacked bar chart and a pie chart. 

> The enterotype analysis is shown in the Enterotype_analysis folder, and a Markdown is also provided to guide the reader step-by-step. The analysis shown is the post-processing analysis after the Encoding/Decoding pipeline (PLEX.I) has ben ran. We refer to Plex.I paper ( https://pubmed.ncbi.nlm.nih.gov/37928248/ ) for details on the Encoding/decoding procedures. 

> "Grouping_LIMMA_results_shared_taxon" folder groups togheter the post-processing for the LIMMA analysi, where the significant taxa identifyied are compared among different cohorts, treatments and time points. 

> ISN_construction folder shown the steps to calculate the ISN starting from a MAGMA network built on population-specific data. 

> LIMMA folder show the LIMMA analysis to find the taxon-taxon interactions that are significantly different between responders and non-responsers. 

> Prediction folder collects the routines and codes to predict responser/non-responder outcomes for each cohort + time + treatment combination.  In detail, the folder is divided into edge-based prediction, where the ISN-edges are used as features, and Matrices_metrics, where are graph metrics, calculated on the edges, that are the predictors.  For both these subfolders, a SVM and a RF routine are implemented. 

> Taxa-level analysis folder only contains the code for the generation of the visual comparison of taxa abundances, for the different phyla, before and after the treatment. 

Finally, scripts Calculate_number_ind_and_taxa_per_cohort.r and Progression_fractional_abundance_test.r contains the code to calculate the number of taxa and individuals for each cohort + time + treatment and the test for the difference of fractional abundance before and after the treatment, respectively. 

As a note, where the repository is provided as "---", the reader needs to substitute with the repository where their data or intermediate results are. 

