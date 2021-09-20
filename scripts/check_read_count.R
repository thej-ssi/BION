#!/usr/bin/env Rscript

library('optparse',warn.conflicts = FALSE)
option_list = list(
  make_option(c("-d", "--otu_table"), type="character", default="otu_table.xlsx",
              help="File name of OTU table in BION format (xlsx or tsv). Default %default", metavar="character"),
  make_option(c("-m", "--metadata_table"), type="character", default="metadata.xlsx",
              help="File name of metadata table in xlsx or tsv format. Default %default", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default="analysis_output",
              help="Name of output directory. Default %default", metavar="character"),
  make_option(c("-v", "--metadata_variable"), type="character", default=NULL,
              help="Name of variable to compare (must match header in metadata table)", metavar="character"),
  make_option(c("-p", "--prokaryot_threshold"), type="integer", default=10000,
              help="Rarefaction threshold for prokaryot reads. Default %default", metavar="integer"),
  make_option(c("-f", "--fungal_threshold"), type="integer", default=1000,
              help="Rarefaction threshold for fungal reads. Default %default", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### Load compact functions ###
print("Loading source functions and packages")
source("/srv/data/tools/git.repositories/BION/scripts/compact_functions_calc.R")

variable_to_compare = opt$metadata_variable # Note that column names in metadata can change when loaded due to illegal characters.

d = load_data(opt$otu_table)
m = load_metadata(opt$metadata_table,metadata_split_variable = variable_to_compare)

dir.create(opt$output_directory, showWarnings = FALSE)

### Check if names in BION output match with names in metadata file ###
compare_datasets(d,m)


### Setup phyloseq objects for Bacterial, eukaryot and fungal data ###

po_all = setup_phylo_object(d,m)

po_pro = po_all$prokaryot
po_eu = po_all$eukaryot
po_fungi = po_all$fungi


test = raw_read_comparison(po_pro,"Geographic.localisation",c())
filename = file.path(opt$output_directory,"raw_read_distribution_prokaryot.png")
ggsave(filename = filename,plot = test$plot,device = "png",dpi = 300, units = "cm", height = 12, width = 20)
