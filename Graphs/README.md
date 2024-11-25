# Graphs
 Collection of all graphical outcomes. 

- In the Alpha_beta diversity subfolder, there are the graphs related to "Alpha Beta diversity" section in the paper. The Alpha_Beta_explanaition.docx document goes into detail. 

- Analysis_Enterotype has the graphs, both in .png and .pdf version, that are whosn in the article for the "Enterotype-based analysis". Additionally, all bact2-bact2 ; bact2-other and other-other graphs are shown. 

- The Enterotype analysis is shown in the Enterotype_analysis folder, and a Markdown is also provided to guide the reader step-by-step. The analysis shown is the post-processing analysis after the Encoding/Decoding pipeline (PLEX.I) has ben ran. We refer to Plex.I paper ( https://pubmed.ncbi.nlm.nih.gov/37928248/ ) for details on the Encoding/decoding procedures. 

- In Microbial_Fractional_abundances there are bar and pie charts showing the abundances of the various cohort + treatment + time, grouped per class, family and phylum. 

- Population_based_connections folder groups the population-based network for each cohort + time + treatment combination, as referred in the "Populaion-based network analysis" of the paper. There networks are both shown based on their taxa, i.e., the lowest level, or the various taxa are grouped into their phyulum and the edge weights averaged. For this, we both took into account the weighted and the binary networks. 

- LIMMA folder shows the LIMMA analysis ("Differential network analysis") to find the taxon-taxon interactions that are significantly different between responders and non-responsers. 

- Prediction folder collects the routines and codes to predict responser/non-responder outcomes for each cohort + time + treatment combination.  In detail, the folder is divided into edge-based prediction, where the ISN-edges are used as features, and Matrices_metrics, where are graph metrics, calculated on the edges, that are the predictors.  For both these subfolders, a SVM and a RF routine are implemented. 

- Taxa-level analysis folder contains the visual comparison of the differences in connectivity (i.e., edge weights) between the same cohort and treatment, but before and after treatment. 
