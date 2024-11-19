

## IT IS not optimal this calling script, but there is no way to have ALL of it scraped
library(stringi)
library(Rcpp)
library(rMAGMA)
library(dplyr)
library(igraph)
library(data.table)

# Parameters --------------------------------------------------------------

eliminate = T # DO we want to eliminate the individuals highlighted in the MAGMA procedure?
diagnosi = c("CD", "UC") ## diagnosi = c("CD","UC")
terapia  =c("VDZ")# c("VDZ") # terapia = c("TNF", "UST", "VDZ")
base_dir = "---"
base_script = "ISN_construction.R"
function_script = "ISN_functions.R"

## ONLY FOR w0
# individual_map = "mapping_nameIndividual_SPARCC.txt"


# OTU_dataset = "w14_n188_OTU_table.txt" # DATASET with the OTU TABLE
OTU_dataset = "w24_n160_OTU_table.txt" # DATASET with the OTU TABLE

# OTU_metadata = "w14_n188_metadata.txt" # METADATA and covariates
OTU_metadata = "w24_n160_metadata.txt" # METADATA and covariates

ASV_annotation = "ASV_taxonomic_annotation.txt"
# weeks = "14"
weeks = "24"

covariates_kept = c("Age", "Gender","FC.nummer", "Disease_duration", "Disease_location_baseline")
identifier_disease = c("FC.nummer","Diagnosis", "Therapy_2" ) # variable in the metadata identifying the cohort and the individuals


setwd(base_dir)
# N.B. BEFORE EVERYTHING IN A NEW FILE, run by hand the first one and save the layout phylum into the global directory 
# e.g. LAYOUT = "./LAYOUT_PHYLUM__graph.txt"
# in that way the graph has (when possible) the same coordinates


eliminate = T; diagnosi = c("CD");terapia  =c("TNF")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)

eliminate = T; diagnosi = c("CD");terapia  =c("UST")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)

eliminate = F; diagnosi = c("CD");terapia  =c("UST")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)

eliminate = T; diagnosi = c("CD");terapia  =c("VDZ")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)


eliminate = T; diagnosi = c("CD", "UC");terapia  =c("TNF")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)

eliminate = T; diagnosi = c("CD", "UC");terapia  =c("VDZ")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)

eliminate = F; diagnosi = c("CD", "UC");terapia  =c("VDZ")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)


eliminate = T; diagnosi = c("UC");terapia  =c("TNF")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)


eliminate = T; diagnosi = c("UC");terapia  =c("VDZ")
setwd(base_dir)
args = list()
args$eliminate = eliminate ; args$terapia = terapia ; args$diagnosi = diagnosi
source (paste0(base_dir, "/Code/ISN_construction/",base_script ))
print(args$eliminate); print(args$diagnosi)
rm(args)
