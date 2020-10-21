### Load compact functions ###
source("https://raw.githubusercontent.com/thej-ssi/BION/master/compact_functions.R")

### Set working directory ###
setwd("/Volumes/data/MPV/THEJ/Projekter/Microbiome/Salem")
setwd("/Users/thej/Documents/GitHub/BION/example/")

### Load BION data (excel file or raw BION tsv output) ###
d = load_data("example_OTU_data.xlsx")

### Load metadata (excel file or tab-separated text file)
m = read.xlsx("example_metadata.xlsx",colNames = TRUE)



### Check if names in BION output match with names in metadata file ###
compare_datasets(d,m)


### Setup phyloseq objects for Bacterial, eukaryot and fungal data ###

po_all = setup_phylo_object(d,m)

po_pro = po_all$prokaryot
po_eu = po_all$eukaryot
po_fungi = po_all$fungi


### Set rarefaction thresholds ###
check_counts(po_pro)
prokaryot_rarefaction_threshold = 11335
rare_pro = rarefy_even_depth(po_pro,prokaryot_rarefaction_threshold,rngseed = 1)

check_counts(po_fungi)
fungal_rarefaction_threshold = 931
rare_fungi = rarefy_even_depth(po_fungi,fungal_rarefaction_threshold,rngseed = 1)



### Indicate which variable in metadata should be compared ###
variable_to_compare = "Geographic.localisation." # Note that column names in metadata can change when loaded due to illegal characters. To see names of metadata variables use:
colnames(sample_data(rare_pro))                  # variable_to_compare should match one of these

color_list = c()  # Colors will be picked automatcially if left empty


### Run prokaryot analysis and save plots ###
name_of_output_folder = "prokaryot_analysis"  # Set name of ouput folder
run_cross_sectional_analysis(rare_pro,variable_to_compare,color_list = color_list,output_folder = name_of_output_folder)

### Run eukaryot analysis and save plots ###
name_of_output_folder = "fungi_analysis"        # Set name of ouput folder
run_cross_sectional_analysis(rare_fungi,variable_to_compare,color_list = color_list,output_folder = name_of_output_folder)
