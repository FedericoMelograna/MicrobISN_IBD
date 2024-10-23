in this folder, 
> Limma_calculation and


Performs the LIMMA, based on ISN, on every cohort/drug/outcome combination. 

After that, Limma_creating_LIMMA_folders_moving_there_all_step2.R moves all the results in a single, common, folder, where, as per the Postprocessing/Computation_significant_pval_LIMMA.R ; we then just filter the results to just take into account hits significant at 0.05 level, after FDR correction. 

