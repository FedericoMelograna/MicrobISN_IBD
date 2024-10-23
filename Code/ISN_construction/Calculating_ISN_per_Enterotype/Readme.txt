In this folder, CALLING_SCRIPT.R calls Enterotype_analysis.R, that perform the enterotype analysis. 
Using routines and scripts previously generated in the general pipeline, we generate the graphs mapping the taxon to the Enterotype. 
In the CALLING_SCRIPT.R, the user need to set the week of analysis ( w0,w14, w24) and the outcome (Clinical, Biomerker, ..). Then, CALLING_SCRIPT.R calls the various combination of cohort/drug. 
Enterotype_analysis.R internally call ISN_construction.R and ISN_functions.R, collectors of functions to build and plot the ISNs. 
