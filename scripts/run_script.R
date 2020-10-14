                ###########################
                #####                 #####
#####################  Set variables  #####################
                #####                 #####
                ###########################

##### BION table and metadata files #####
setwd("P:/MPV/THEJ/Projekter/Microbiome/Biller_201802/") ### NEEDED INPUT. Folder containing input files
input_file = "all_data_nobeetle.xlsx"                    ### NEEDED INPUT. Name of file containing sequence counts (.txt og .xlsx)
metadata_file = "metadata.txt"                           ### NEEDED INPUT. Name of metadata file (.txt og .xlsx)

##### metadata variable for sub set plots #####
metadata_split_variable = "type"                         ### NEEDED INPUT. Header of the column that indicates which group the sample is in ###
metadata_split_groups = c("I2","I7","U2","U7")           ### NEEDED INPUT. Names of the different groups to be compared



### Load source file ###
#source("P:/MPV/THEJ/BION/FABIAN_backup/automate_functions.R")
source("https://raw.githubusercontent.com/thej-ssi/BION/master/automate_functions.R")

##### Load input files #####
tsv_input = load_data(input_file)
metadata = load_metadata(metadata_file,metadata_split_variable)

### Create phyloseq objects for prokaryot, species, eukaryot species, and funal species ###
phylo_objects<-setup_phylo_object_2(tsv_input,metadata)
po_prokaryot<-phylo_objects[[1]]
po_eukaryot<-phylo_objects[[2]]
po_fungi<-phylo_objects[[3]]

                                                      
############################################################################
################### Find suitable rarefaction thresholds ###################
############################################################################

check_counts(po_prokaryot)
check_counts(po_eukaryot)
check_counts(po_fungi)


##### Set rarefaction thresholds #####
prokaryot_rarefaction_threshold = 26040         ### NEEDED INPUT, rarefaction threshold for prokaryotes
eukaryot_rarefaction_threshold = 5760           ### NEEDED INPUT, rarefaction threshold for eukaryotes
fungi_rarefaction_threshold = 1560              ### NEEDED INPUT, rarefaction threshold for fungi


################################
#####                      #####
#####    Function calls    #####
#####                      #####
################################


### Rarefy ###
rare_prokaryot = rarefy_even_depth(po_prokaryot,sample.size = prokaryot_rarefaction_threshold,rngseed=1)
rare_eukaryot = rarefy_even_depth(po_eukaryot,sample.size = eukaryot_rarefaction_threshold,rngseed=1)
rare_fungi = rarefy_even_depth(po_fungi,sample.size = fungi_rarefaction_threshold,rngseed=1)


### Run all analyses ###
summary_list = list(prokaryot_rarefaction_threshold,eukaryot_rarefaction_threshold,fungi_rarefaction_threshold,nsamples(po_prokaryot),nrow(tsv_input))
run_all(rare_prokaryot,rare_eukaryot,rare_fungi,summary_list)


