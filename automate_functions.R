
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
  require(x,character.only = TRUE)
}


##### Check if dependencies are installed and load packages #####
if(! require("RCurl")) {install.packages("RCurl")}
library('RCurl')

if(! require("phyloseq")) {
  source('http://bioconductor.org/biocLite.R')
  biocLite('phyloseq')
}
library('phyloseq')

if(! require("DESeq2")) {
  source('http://bioconductor.org/biocLite.R')
  biocLite("DESeq2")
}
library('DESeq2')

pkgTest('vegan')
pkgTest('tibble')
pkgTest('reshape')
pkgTest('plotly')
pkgTest('gplots')
pkgTest('yaml')
pkgTest('backports')
pkgTest('openxlsx')
pkgTest('webshot')
pkgTest('stringr')
pkgTest('plyr')

#webshot::install_phantomjs()


setup_outdir <- function() {
  n = 1
  flag = 0
  while (flag == 0) {
    outdir = paste0('FABIAN_output_',n)
    if (file.exists(outdir)) {
      n = n+1
    } else {
      flag = 1
      dir.create(outdir)
    }
  }
  return(outdir)
}




load_data <- function(input_file) {
  if (tolower(substr(input_file,nchar(input_file)-4,nchar(input_file))) == ".xlsx") {
    tsv_input = read.xlsx(input_file, 1) # read the first sheet
  }
  else {
    tsv_input = read.table(input_file,sep = '\t', comment.char = "", header = TRUE, stringsAsFactors = FALSE,check.names = FALSE)
    if (substr(tsv_input[1,1],1,1) == '#') {
      colcount = ncol(tsv_input)
      sample_vector = colnames(tsv_input)[c(1:(colcount-6))]
      sample_vector[1] = substr(sample_vector[1],3,nchar(sample_vector[1]))
      counts = tsv_input[c(4:nrow(tsv_input)),c(1:(colcount-6))]
      tax_vector = as.vector(tsv_input[c(4:nrow(tsv_input)),ncol(tsv_input)])
      colnames(counts) = sample_vector
      
      new_counts = matrix(nrow = nrow(counts), ncol = 0)
      for (col in 1:ncol(counts)) {
        new_vector = as.numeric(as.vector(counts[,col]))
        new_counts = cbind(new_counts,new_vector)
      }
      colnames(new_counts) = sample_vector
      tsv_input = as.data.frame(new_counts)
      tsv_input$Taxonomic.groups = tax_vector
    }
  }
  return(tsv_input)
}

BIONtsv_2_phylodata <- function(input_file) {
  if (tolower(substr(input_file,nchar(input_file)-4,nchar(input_file))) == ".xlsx") {
    tsv_input = read.xlsx(input_file, 1) # read the first sheet    tax_vector = as.vector(tsv_input$Taxonomic.groups)
    
  }
  else {
    tsv_input = read.table(input_file,sep = '\t', comment.char = "", header = TRUE, stringsAsFactors = FALSE,check.names = FALSE)
    if (substr(tsv_input[1,1],1,1) == '#') {
      colcount = ncol(tsv_input)
      sample_vector = colnames(tsv_input)[c(1:(colcount-6))]
      sample_vector[1] = substr(sample_vector[1],3,nchar(sample_vector[1]))
      counts = tsv_input[c(4:nrow(tsv_input)),c(1:(colcount-6))]
      tax_vector = as.vector(tsv_input[c(4:nrow(tsv_input)),ncol(tsv_input)])
      colnames(counts) = sample_vector
      
      new_counts = matrix(nrow = nrow(counts), ncol = 0)
      for (col in 1:ncol(counts)) {
        new_vector = as.numeric(as.vector(counts[,col]))
        new_counts = cbind(new_counts,new_vector)
      }
      colnames(new_counts) = sample_vector
      tsv_input = as.data.frame(new_counts)
      tsv_input$Taxonomic.groups = tax_vector
    }
    otu_table = tsv_input[,!colnames(tsv_input) %in% c("Taxonomic.groups","Row.max")]
    tax_table = matrix(nrow = 0, ncol = 7)
    for (i in 1:length(tax_vector)) {
      tax_split = strsplit(tax_vector[i],"; ",fixed=TRUE)
      tax_vec = c(tax_split[[1]][1],tax_split[[1]][2],tax_split[[1]][3],tax_split[[1]][4],tax_split[[1]][5],tax_split[[1]][6],tax_split[[1]][7])
      tax_vec = unlist(lapply(tax_vec,function(x) substr(x,4,nchar(x))))
      tax_table = rbind(tax_table,tax_vec)
    }
  }
  colnames(tax_table) = c("Domain","Phylum","Class","Order","Family","Genus","Species")
  rownames(otu_table) = paste0("OTU",1:nrow(otu_table))
  rownames(tax_table) = rownames(otu_table)
  return(list(tax.table=tax_table,otu.table=otu_table))
}



load_metadata <- function(metadata_file, metadata_split_variable) {
  if (metadata_file == "" | is.na(metadata_file)) {
    metadata = ""
  } else {
    if (tolower(substr(metadata_file,nchar(metadata_file)-4,nchar(metadata_file))) == ".xlsx") {
      metadata = read.xlsx(metadata_file, 1,rowNames = TRUE) # read the first sheet
    }
    else {
      metadata = read.table(metadata_file,sep = '\t', header = TRUE, row.names = 1)
    }
    sample_groups = as.character(metadata[,colnames(metadata)==metadata_split_variable])
    metadata = cbind(metadata,sample_groups)
    rownames(metadata) = gsub("-",".",rownames(metadata),)
    #metadata$sample_groups = gsub("-",".",metadata$sample_groups)
  }
  return(metadata)
}



collapse_Aspergillaceae <- function(tax_table,otu_table) {
  Aspergillaceae_rows = which(tax_table[,5]=="Aspergillaceae")
  Trichocomaceae_rows = which(tax_table[,5]=="Trichocomaceae")
  for (i in Aspergillaceae_rows) {
    asper_tax = as.vector(tax_table[i,])
    asper_otu = as.vector(otu_table[i,])
    asper_genus = asper_tax[6]
    asper_species = asper_tax[7]
    remove_rows = c()
    for (j in Trichocomaceae_rows) {
      tricho_tax = as.vector(tax_table[j,])
      tricho_otu = as.vector(otu_table[j,])
      tricho_genus = tricho_tax[6]
      tricho_species = tricho_tax[7]
      matchflag = 0
      if (asper_genus==tricho_genus & asper_species==tricho_species) {
        new_otu = asper_otu+tricho_otu
        otu_table[j,] = new_otu
        remove_rows = c(remove_rows,i)
        matchflag = 1
      }
      if (matchflag == 0) {
        tax_table[i,]=c(as.vector(tax_table[i,1:4]),"Trichocomaceae",as.vector(tax_table[i,6:7]))
      }
      
    }
    if (length(remove_rows)>0){
      tax_table = tax_table[-remove_rows,]
      otu_table = otu_table[-remove_rows,]
    }
  }
  Aspergillaceae_rows = which(tax_table[,5]=="Aspergillaceae")
  for (i in Aspergillaceae_rows) {
    asper_tax = as.vector(tax_table[i,])
    asper_tax[5] = 'Trichocomaceae'
    tax_table[i,] = asper_tax
    print(asper_tax)
  }
  return_list = list(tax_table,otu_table)
}

remove_duplicates_from_phylo_object_2 <- function(po) {
  tax_string<-tsv_input[,1] # din out tabel, med taxonomy i første kolonne
  tax_string<-gsub("\\(","",tax_string)
  tax_string<-gsub("\\)","",tax_string)
  tsv_input[,1]<-tax_string
  tsv_input_summed<-ddply(tsv_input, .(Taxonomic.groups), numcolwise(sum))
}




combine_duplicates_from_phylo_object <- function(po) {
  d = otu_table(po)
  tax = tax_table(po)
  new_d = matrix(nrow=0,ncol = ncol(d))
  new_tax = matrix(nrow=0,ncol = ncol(tax))
  otu_vec = c()
  all_tax_strings = c()
  for (row in 1:nrow(tax)) {
    tax_vec = c()
    for (col in 1:ncol(tax)) {
      old_tax_pos = as.character(tax[row,col])
      if (grepl('%',old_tax_pos)) {
        tax_vec = c(tax_vec,'NA')
      } else if (grepl('\\(',old_tax_pos)) {
        new_tax_pos = substr(old_tax_pos,2,(nchar(old_tax_pos)-1))
        #print(new_tax_pos)
        tax_vec = c(tax_vec,new_tax_pos)
        #print(tax_vec)
        
      } else {
        tax_vec = c(tax_vec,tax[row,col])
      }
    }
    tax_string = paste(tax_vec, collapse = '|')
    pos = match(tax_string, all_tax_strings)
    if (!is.na(pos)) {
      new_d[pos,] = new_d[pos,]+d[row,]
    } else {
      new_d = rbind(new_d,d[row,])
      new_tax = rbind(new_tax,tax_vec)
      otu_vec = c(otu_vec,rownames(tax)[row])
      all_tax_strings = c(all_tax_strings,tax_string)
    }
  }
  rownames(new_d) = otu_vec
  rownames(new_tax) = otu_vec
  colnames(new_tax) = colnames(tax)
  return_po = phyloseq(tax_table(new_tax),otu_table(new_d,taxa_are_rows = TRUE),sample_data(sample_data(po)))
  return(return_po)
}



#######################################################
##### Merge all raw BION output probe files       #####
##### Set all counts less than x% of row max to 0 #####
##### Combine all Mammalia                        #####
##### Remove bacteria from 18s probe files        #####
##### Collapse suggested taxa to unclassified     #####
#######################################################

denoise_microbiome <- function(data_tables, threshold, data_type) {
  ## read in files, and sort for librarysize and count
  final<-list()
  i<-1
  for(data_table in data_tables){
    
    sample_data<-data_table #read the sample
    #sample_data<-data_tables[[1]] #REMOVE, only for testing
    type<-data_type[i]
    
    # remove unecessary columns
    columns_to_remove<-c("Row.max","Row.sum", "Sim.", "Fav.", "Mark")
    sample_data<-sample_data[,!(colnames(sample_data) %in% columns_to_remove)] # remove row max etc. 
    
    # replace >=threshold with 0, and rows that thereafter equals zero
    sample_data[sample_data<=threshold]<-0 # set everything below or equal to 3 to zero, do his after merging instead? Should not matter if I merge with max inseat of sum...
    sample_data$new_row_max<-apply(sample_data[1:(ncol(sample_data)-1)], 1, max) # remove rows where row_max=0
    sample_data_max_over_zero<-sample_data[sample_data$new_row_max>0,]
    sample_data<-sample_data_max_over_zero[,-ncol(sample_data_max_over_zero)] # remove sum column
    names(sample_data)[ncol(sample_data)]<-"Taxonomic.groups" # change names of the last column to "Taxonomic.groups" (skipping line numbers) 
    
    #remove bacteria from 18S data
    print(type)
    print(paste("number of OTUs before removing bacteria:",toString(nrow(sample_data))))
    if(type=="18S"){
      rem_lines_index<-grep("d__Bacteria;", sample_data$Taxonomic.groups)
      if(length(rem_lines_index>0)){
        sample_data<-sample_data[-rem_lines_index,]
      }
    }
    print(paste("number of OTUs after removing bacteria:",toString(nrow(sample_data))))
    
    ### should merge twice, ones to merge primer sets and one to merge after non-favorites are moved up one level. Make correct names sort with sum per dataset, THEN merge primersets with max
    
    ### add "unclassified" if values are missing, e.g. stopping at genus instead of species
    tax_strings<-sample_data$Taxonomic.groups
    count <- str_count(tax_strings,";")
    classifier=c("; p__unclassified","; c__unclassified","; o__unclassified","; f__unclassified","; g__unclassified","; s__unclassified")
    for(ii in 0:length(classifier)-1){
      tax_strings[count==ii] = paste0(tax_strings[count==ii],classifier[ii+1])
      count <- str_count(tax_strings,";")
    }
    
    tax_frame<-colsplit(tax_strings, "; ", names=c("domain","phylum","class","order","family","genus","species"))
    
    ### go one up in hiearcy in cases of non-favorites
    taxmat_d<-data.frame("kingdom"=as.character(tax_frame[,1]), "phylum"=as.character(tax_frame[,2]), "class"=as.character(tax_frame[,3]), "order"=as.character(tax_frame[,4]), "family"=as.character(tax_frame[,5]), "genus"=as.character(tax_frame[,6]), "species"=as.character(tax_frame[,7]), stringsAsFactors=F)
    
    taxmat_d$species[grep("%",taxmat_d$species)]<-"s__unclassified"
    taxmat_d$genus[grep("%",taxmat_d$genus)]<-"g__unclassified"
    taxmat_d$family[grep("%",taxmat_d$family)]<-"f__unclassified"
    taxmat_d$order[grep("%",taxmat_d$order)]<-"o__unclassified"
    taxmat_d$class[grep("%",taxmat_d$class)]<-"c__unclassified"
    taxmat_d$phylum[grep("%",taxmat_d$phylum)]<-"p__unclassified"
    taxmat_d$kingdom[grep("%",taxmat_d$kingdom)]<-"d__unclassified"
    new_tax_names<-apply(taxmat_d, 1, paste, collapse="; ")
    
    sample_data$Taxonomic.groups<-new_tax_names
    
    # Merge within primer set by summing (that is: e.g. unclassified and non-favorites are now called the same e.g.  staph_unclassified, with staph_aureus60%/epi30%/lug10&)
    sample_data_summed<-ddply(sample_data, .(Taxonomic.groups), numcolwise(sum))
    #print(paste("number of rows before merging:",toString(nrow(sample_data))))
    print(paste("number of rows after merging:",toString(nrow(sample_data_summed))))
    
    final<-rbind(final, sample_data_summed)
    i<-i+1
  }
  ### remove mammalia
  #mammal_rem_index<-grep("Mammalia", taxmat_d$class)
  #print(paste("number of mammalian sequences removed:",toString(length(mammal_rem_index))))
  #final_mammals_rem<-final[-mammal_rem_index,]
  
  print("## further sorting ##")
  
  ### Merge mammals
  tax_frame<-colsplit(final$Taxonomic.groups, "; ", names=c("domain","phylum","class","order","family","genus","species"))
  mammal_index<-grep("Mammalia", tax_frame$class)
  if(length(mammal_index)>0){
    mammals_table<-final[mammal_index,]
    mammals_identical_merged<-ddply(mammals_table, .(Taxonomic.groups), numcolwise(max))  # several identical lines, merge these by tkaing max (contribution from diff. primer set)
    mammals_summed<-apply(mammals_identical_merged[2:ncol(mammals_identical_merged)], 2, sum) # sum up the rest (Tax_group, shifts to column one using ddply)
    mammals_summed<-as.data.frame(t(mammals_summed))
    mammals_summed$Taxonomic.groups<-"d__Eukaryota; p__Vertebrata; c__Mammalia; o__unclassified; f__unclassified; g__unclassified; s__unclassified"
    mammals_summed<-mammals_summed[,c(ncol(mammals_summed), 1:(ncol(mammals_summed)-1))]
    final_mammals_rem<-final[-mammal_index,]  # Remove mammals from the original final frame, and add the new mammal line (name of dataframe misleading...)
    final_mammals_rem<-rbind(final_mammals_rem, mammals_summed)
    print(paste("all mammalian: ", toString(nrow(mammals_table))))
    print(paste("number of mammalian sequences that have been merged: ", toString(nrow(mammals_identical_merged))))
  }else{
    final_mammals_rem<-final
  }
  
  ### merge identical OTUs
  print(paste("number of rows before mering:",toString(nrow(final_mammals_rem))))
  final_summed<-ddply(final_mammals_rem, .(Taxonomic.groups), numcolwise(max))
  print(paste("number of rows after mering:",toString(nrow(final_summed))))
  
  return(final_summed)
  
}

data_initialization <- function(input_file_list) { ## Add this to Thors code
  files_found=TRUE
  
  for(input_file in input_file_list) {
    if(!file.exists(file.path(input_file))) {
      cat(sprintf("%s not found\n",input_file))
      files_found=FALSE
    }
  }
  if(!files_found) {
    stop("Missing files")
  }
  data_tables <- list()
  for(input_file in input_file_list){
    data <- readLines(input_file)
    if(grep(("#\\s*"),data[1])==1){ # if first line start with hash tag, remove hash tag
      data[1] <- gsub("#\\s*","",data[1])
    }
    else{
      stop(cat(sprintf("%s is not in proper format\n",input_file)))
    }
    data_tables<-c(data_tables,list(read.table(text=data, sep="\t", h=T, stringsAsFactors=F)))
  }
  
  print(paste("number of data tables added:",toString(length(data_tables))))
  return(data_tables)
  
}




### Set up two phyloseq objects, one for prokaryot and one for eukaryot species ###
setup_phylo_object <- function(tsv_input,metadata=NA) {
  if (colnames(tsv_input)[1] == "domain") {
    otu_table = as.matrix(tsv_input[,-(1:7)])
    tax_table = as.matrix(tsv_input[,1:7])
  } else {
    tax_vector = as.vector(tsv_input$Taxonomic.groups)
    otu_table = tsv_input[,!colnames(tsv_input) %in% c("Taxonomic.groups","Row.max")]
    tax_table = matrix(nrow = 0, ncol = 7)
    for (i in 1:length(tax_vector)) {
      tax_split = strsplit(tax_vector[i],"; ",fixed=TRUE)
      tax_vec = c(tax_split[[1]][1],tax_split[[1]][2],tax_split[[1]][3],tax_split[[1]][4],tax_split[[1]][5],tax_split[[1]][6],tax_split[[1]][7])
      tax_table = rbind(tax_table,tax_vec)
    }
  }
  tax_table = trimws(tax_table)
  tax_table = substr(tax_table,4,nchar(tax_table))
  Aspergillaceae_remove = collapse_Aspergillaceae(tax_table,otu_table)
  tax_table = Aspergillaceae_remove[[1]]
  otu_table = Aspergillaceae_remove[[2]]
  tax_prokaryot = as.matrix(tax_table[(tax_table[,1]=="Bacteria"|tax_table[,1]=="Archaea") & !tax_table[,2]=="Cyanobacteria/Chloroplast",])
  otu_prokaryot = data.matrix(otu_table[(tax_table[,1]=="Bacteria"|tax_table[,1]=="Archaea") & !tax_table[,2]=="Cyanobacteria/Chloroplast",])
  tax_eukaryot = as.matrix(tax_table[tax_table[,1]=="Eukaryota" & !tax_table[,3]=="Mammalia",])
  otu_eukaryot = data.matrix(otu_table[tax_table[,1]=="Eukaryota" & !tax_table[,3]=="Mammalia",])
  fungi_rows = c()
  for (n in 1:nrow(tax_table)) {
    phylum = as.character(tax_table[n,2])
    if (!is.na(phylum) & nchar(phylum)>=8) {
      #print(phylum)
      substr1 = substr(phylum,nchar(phylum)-5,nchar(phylum))
      substr2 = substr(phylum,nchar(phylum)-7,nchar(phylum))
      if(substr1=="mycota" | substr2=="mycotina") {
        fungi_rows = c(fungi_rows,n)
      }
    }
  }
  tax_ranks = c("Domain","Phylum","Class","Order","Family","Genus","Species")
  if (metadata == "" | is.na(metadata)) {
    metadata = cbind(rep("all samples"),ncol(otu_prokaryot))
    rownames(metadata) = colnames(otu_prokaryot)
    colnames(metadata)[1] = "sample_groups"
  }
  if (nrow(tax_prokaryot) > 0) {
    rownames(tax_prokaryot) <- paste0("OTU", 1:nrow(tax_prokaryot))
    rownames(otu_prokaryot) = rownames(tax_prokaryot)
    colnames(tax_prokaryot) = tax_ranks
    po_prokaryot = phyloseq(tax_table(tax_prokaryot),otu_table(otu_prokaryot,taxa_are_rows = TRUE),sample_data(metadata))
  } else {
    po_prokaryot = 0
  }
  if (nrow(tax_eukaryot) > 0) {
    rownames(tax_eukaryot) <- paste0("OTU", 1:nrow(tax_eukaryot))
    rownames(otu_eukaryot) = rownames(tax_eukaryot)
    colnames(tax_eukaryot) = tax_ranks
    po_eukaryot = phyloseq(tax_table(tax_eukaryot),otu_table(otu_eukaryot,taxa_are_rows = TRUE),sample_data(metadata))
  } else {
    po_eukaryot = 0
  }
  if (length(fungi_rows)>0) {
    tax_fungi = tax_table[fungi_rows,]
    otu_fungi = otu_table[fungi_rows,]
    rownames(tax_fungi) <- paste0("OTU", 1:nrow(tax_fungi))
    rownames(otu_fungi) = rownames(tax_fungi)
    colnames(tax_fungi) = tax_ranks
    po_fungi = phyloseq(tax_table(tax_fungi),otu_table(otu_fungi,taxa_are_rows = TRUE),sample_data(metadata))
  } else {
    po_fungi = 0
  }
  returnlist = list('prokaryot'=po_prokaryot,'eukaryot'=po_eukaryot,'fungi'=po_fungi)
  return(returnlist)
}


### Set up two phyloseq objects, one for prokaryot and one for eukaryot species ###
setup_phylo_object_2 <- function(tsv_input,metadata) {
  otu_table = as.matrix(tsv_input[,-(1:7)])
  tax_table = as.matrix(tsv_input[,1:7])
  tax_table = trimws(tax_table)
  tax_table = substr(tax_table(po),4,nchar(tax_table(po)))
  tax_prokaryot = tax_table[tax_table[,1]=="Bacteria"|tax_table[,1]=="Archaea",]
  otu_prokaryot = otu_table[tax_table[,1]=="Bacteria"|tax_table[,1]=="Archaea",]
  tax_eukaryot = tax_table[tax_table[,1]=="Eukaryota" & !tax_table[,3]=="Mammalia",]
  otu_eukaryot = otu_table[tax_table[,1]=="Eukaryota" & !tax_table[,3]=="Mammalia",]
  tax_ranks = c("Domain","Phylum","Class","Order","Family","Genus","Species")
  rownames(tax_prokaryot) <- paste0("OTU", 1:nrow(tax_prokaryot))
  rownames(otu_prokaryot) = rownames(tax_prokaryot)
  rownames(tax_eukaryot) <- paste0("OTU", 1:nrow(tax_eukaryot))
  rownames(otu_eukaryot) = rownames(tax_eukaryot)
  colnames(tax_prokaryot) = tax_ranks
  colnames(tax_eukaryot) = tax_ranks
  po_prokaryot = phyloseq(tax_table(tax_prokaryot),otu_table(otu_prokaryot,taxa_are_rows = TRUE),sample_data(metadata))
  po_eukaryot = phyloseq(tax_table(tax_eukaryot),otu_table(otu_eukaryot,taxa_are_rows = TRUE),sample_data(metadata))
  returnlist = list(po_prokaryot,po_eukaryot)
  return(returnlist)
}


### Create individual phyloseq for subsets of samples specified by metadata ###
split_by_metadata <- function(po,metadata_split_variable,metadata_split_groups) {
  variable_vector = as.vector(sample_data(po)[,colnames(sample_data(po))==metadata_split_variable])
  return_po = list()
  for (i in 1:length(metadata_split_groups)) {
    samples = rownames(sample_data(po))[variable_vector==metadata_split_groups[i]]
    if (length(samples)>0) {
      temp_po = prune_samples(samples,po)
    } else {
      temp_po = 0
    }
    return_po = c(return_po,temp_po)
  }
  return(return_po)
}

get_top10_taxa <- function(po) {
  otu_sum<-apply(otu_table(po), 1, sum)
  otus_sorted<-sort.list(otu_sum, decreasing=T)
  top10_otu=otus_sorted[1:10]
  return(top10_otu)
}

get_top_n_taxa <- function(po,n) {
  otu_sum<-apply(otu_table(po), 1, sum)
  otus_sorted<-sort.list(otu_sum, decreasing=T)
  top_n_otu=otus_sorted[1:n]
  return(top_n_otu)
}

get_taxa_names <- function(po,taxa) {
  taxmat = tax_table(po)
  newnames = c()
  for (i in 1:length(taxa)) {
    rownumber = taxa[i]
    tax_vector = as.vector(taxmat[rownumber,])
    if (!is.na(tax_vector[7])) {
      newname = paste0(tax_vector[6],' ',tax_vector[7])
    } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
      newname = tax_vector[6]
    } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
      newname = tax_vector[5]
    } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
      newname = tax_vector[4]
    } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
      newname = tax_vector[3]
    } else {
      newname = tax_vector[2]
    }
    newnames = c(newnames,newname)
  }
  return(newnames)
}

make_abundance_barplot <- function(po,taxa,plot_name="Relative abundance") {
  if (class(po)=="phyloseq") {
    taxmat = tax_table(po)
    dd<-otu_table(po)
    dd<-apply(dd, 2, function(x) x/sum(x)*100)
    dd<-as.data.frame(dd)
    
    dd$sum<-apply(dd, 1, sum)
    dd_sorted<-dd[dd$sum>0,]
    dd_sorted<-dd_sorted[,-ncol(dd_sorted)]
    dd_sorted_d<-vegdist(t(dd_sorted), method="bray")
    fit <- hclust(dd_sorted_d, method="ward.D")
    #plot(fit, cex=0.5) # display dendogram
    cluster_order<-fit$labels[fit$order]
    newnames = c()
    top10_otu_table = matrix(nrow=0,ncol=ncol(dd_sorted))
    for (i in 1:length(taxa)) {
      rownumber = taxa[i]
      tax_vector = as.vector(taxmat[rownumber,])
      if (!is.na(tax_vector[7])) {
        newname = paste0(tax_vector[6],' ',tax_vector[7])
      } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
        newname = tax_vector[6]
      } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
        newname = tax_vector[5]
      } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
        newname = tax_vector[4]
      } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
        newname = tax_vector[3]
      } else {
        newname = tax_vector[2]
      }
      newnames = c(newnames,newname)
      top10_otu_table = rbind(top10_otu_table,dd[rownumber,])
    }
    rownames(top10_otu_table)<-newnames
    top10=top10_otu_table
    top10 = top10[,!colnames(top10)=="sum"]
    top10$species<-row.names(top10)
    melt_top10<-melt(top10)
    melt_top10$variable <- factor(melt_top10$variable,levels = cluster_order)
    names(melt_top10)<-c("genus", "ID", "percent")
    p <- plot_ly(data = melt_top10[which(melt_top10$genus == rownames(top10)[10]),], x = ~ID, y = ~percent, type = 'bar', name =rownames(top10)[10]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[9]),], name =rownames(top10)[9]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[8]),], name =rownames(top10)[8]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[7]),], name =rownames(top10)[7]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[6]),], name =rownames(top10)[6]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[5]),], name =rownames(top10)[5]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[4]),], name =rownames(top10)[4]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[3]),], name =rownames(top10)[3]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[2]),], name =rownames(top10)[2]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[1]),], name =rownames(top10)[1]) %>%
      layout(title = plot_name,
             xaxis=list(title="Sample ID"),
             yaxis=list(title="Abundance (percent of rarefied counts)"),
             barmode = 'stack',
             #autosize = F,
             margin = list(l=50,r=50,b=100,t=50))
    p
    return(p)
  }
}

make_abundance_barplot_2 <- function(po,top10_species,plot_name) {
  if (class(po)=="phyloseq") {
    taxmat = tax_table(po)
    dd<-otu_table(po)
    dd<-apply(dd, 2, function(x) x/sum(x)*100)
    dd<-as.data.frame(dd)
    
    dd$sum<-apply(dd, 1, sum)
    dd_sorted<-dd[dd$sum>0,]
    dd_sorted<-dd_sorted[,-ncol(dd_sorted)]
    dd_sorted_d<-vegdist(t(dd_sorted), method="bray")
    fit <- hclust(dd_sorted_d, method="ward.D")
    #plot(fit, cex=0.5) # display dendogram
    cluster_order<-fit$labels[fit$order]
    newnames = c()
    top10_otu_table = matrix(nrow=0,ncol=ncol(dd_sorted))
    for (i in 1:length(top10_species)) {
      rownumber = top10_species[i]
      tax_vector = as.vector(taxmat[rownumber,])
      if (!is.na(tax_vector[7])) {
        newname = paste0(tax_vector[6],' ',tax_vector[7])
      } else if (!tax_vector[6]=="unclassified") {
        newname = tax_vector[6]
      } else if (!tax_vector[5]=="unclassified") {
        newname = tax_vector[5]
      } else if (!tax_vector[4]=="unclassified") {
        newname = tax_vector[4]
      } else if (!tax_vector[3]=="unclassified") {
        newname = tax_vector[3]
      } else {
        newname = tax_vector[2]
      }
      newnames = c(newnames,newname)
      top10_otu_table = rbind(top10_otu_table,dd[rownumber,])
    }
    rownames(top10_otu_table)<-newnames
    top10=top10_otu_table
    top10 = top10[,!colnames(top10)=="sum"]
    top10$species<-row.names(top10)
    melt_top10<-melt(top10)
    melt_top10$variable <- factor(melt_top10$variable,levels = cluster_order)
    names(melt_top10)<-c("genus", "ID", "percent")
    p <- plot_ly(data = melt_top10[which(melt_top10$genus == rownames(top10)[10]),], x = ~ID, y = ~percent, type = 'bar', name =rownames(top10)[10]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[9]),], name =rownames(top10)[9]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[8]),], name =rownames(top10)[8]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[7]),], name =rownames(top10)[7]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[6]),], name =rownames(top10)[6]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[5]),], name =rownames(top10)[5]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[4]),], name =rownames(top10)[4]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[3]),], name =rownames(top10)[3]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[2]),], name =rownames(top10)[2]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[1]),], name =rownames(top10)[1]) %>%
      layout(title = plot_name,
             xaxis=list(title="Sample ID"),
             yaxis=list(title="Abundance (percent of rarefied counts)"),
             barmode = 'stack',
             #autosize = F,
             margin = list(l=50,r=50,b=100,t=50))
    p
    return(p)
  }
}

alphadiversity_plot <- function(po,plot_name) {
  col_vec = RColorBrewer::brewer.pal(length(levels(factor(get_variable(po)$sample_groups))),"Set1")
  r <- data.frame(ID=sample_names(po), type=factor(get_variable(po)$sample_groups), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon")))
  p1 <- plot_ly(r, y = ~Observed, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5, colors = col_vec) %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = 90),
           yaxis=list(title='Number of observed OTUs'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  p2 <- plot_ly(r, y = ~Shannon, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5,colors = col_vec) %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = -45),
           yaxis=list(title='Shannon diversity index'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  kruskal_Observed = kruskal.test(r$Observed,r$type)
  kruskal_Shannon = kruskal.test(r$Shannon,r$type)
  print(p1)
  print(p2)
  return_list = list("Observed_plot"=p1,"Shannon_plot"=p2,
                     "Observed_pvalue"=kruskal_Observed$p.value,"Shannon_pvalue"=kruskal_Shannon$p.value)
  return(return_list)
}

make_alphadiversity_object <- function(po,variable_name,plot_title,color_list) {
  groups = levels(factor(get_variable(po,variable_name)))
  if (length(color_list) == length(groups)) {
    col_vec = color_list
  } else {
    col_vec = RColorBrewer::brewer.pal(length(groups),"Set1")
    print(paste0('Number of colors given (', length(color_list) , ') does not match number of levels in variable (', length(groups),')'))
  }
  r <- data.frame(ID=sample_names(po), type=factor(get_variable(po,variable_name)), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon")))
  p1 <- plot_ly(r, y = ~Observed, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5, colors = col_vec) %>%
    layout(title = plot_title,
           #xaxis=list(tickangle = 90),
           yaxis=list(title='Number of observed OTUs'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  p2 <- plot_ly(r, y = ~Shannon, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5,colors = col_vec) %>%
    layout(title = plot_title,
           #xaxis=list(tickangle = -45),
           yaxis=list(title='Shannon diversity index'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  shannon_matrix = matrix(ncol = length(groups), nrow = length(groups))
  observed_matrix = matrix(ncol = length(groups), nrow = length(groups))
  for (n1 in 1:(length(groups)-1)) {
    print(n1)
    for (n2 in (n1+1):length(groups)) {
      group1 = groups[n1]
      group2 = groups[n2]
      ### Shannon
      print(paste0('group1 ',group1))
      print(paste0('group2 ',group2))
      vec1 = r$Shannon[which(r$type==group1)]
      vec2 = r$Shannon[which(r$type==group2)]
      wilcox_shannon = wilcox.test(vec1,vec2)
      shannon_matrix[n1,n2] <- wilcox_shannon$p.value
      shannon_matrix[n2,n1] <- wilcox_shannon$p.value
      ### Observed
      print(paste0('group1 ',group1))
      print(paste0('group2 ',group2))
      vec1 = r$Observed[which(r$type==group1)]
      vec2 = r$Observed[which(r$type==group2)]
      wilcox_observed = wilcox.test(vec1,vec2)
      observed_matrix[n1,n2] <- wilcox_observed$p.value
      observed_matrix[n2,n1] <- wilcox_observed$p.value
    }
    
  }
  kruskal_Observed = kruskal.test(r$Observed,r$type)
  kruskal_Shannon = kruskal.test(r$Shannon,r$type)
  p_df_shannon = as.data.frame(shannon_matrix)
  p_df_observed = as.data.frame(observed_matrix)
  rownames(p_df_shannon) = groups
  colnames(p_df_shannon) = groups
  rownames(p_df_observed) = groups
  colnames(p_df_observed) = groups
  return_list = list("Observed_plot"=p1,"Shannon_plot"=p2,
                     "Observed_kruskal"=kruskal_Observed$p.value,"Shannon_kruskal"=kruskal_Shannon$p.value,
                     "Observed_MWU_mat"=p_df_observed,"Shannon_MWU_mat"=p_df_shannon)
  print(col_vec)
  return(return_list)
}
make_barplot_with_tiles <- function(po,taxa,variable_name,plot_title="",color_vector=c()) {
  if (class(po)=="phyloseq") {
    taxmat = tax_table(po)
    dd<-otu_table(po)
    dd<-apply(dd, 2, function(x) x/sum(x)*100)
    dd<-as.data.frame(dd)
    split_variable_vector = get_variable(po,variable_name)
    dd$sum<-apply(dd, 1, sum)
    dd_sorted<-dd[dd$sum>0,]
    dd_sorted<-dd_sorted[,-ncol(dd_sorted)]
    dd_sorted_d<-vegdist(t(dd_sorted), method="bray")
    fit <- hclust(dd_sorted_d, method="ward.D")
    #plot(fit, cex=0.5) # display dendogram
    cluster_order<-fit$labels[fit$order]
    newnames = c()
    top10_otu_table = matrix(nrow=0,ncol=ncol(dd_sorted))
    for (i in 1:length(taxa)) {
      rownumber = taxa[i]
      tax_vector = as.vector(taxmat[rownumber,])
      if (!is.na(tax_vector[7])) {
        newname = paste0(tax_vector[6],' ',tax_vector[7])
      } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
        newname = tax_vector[6]
      } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
        newname = tax_vector[5]
      } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
        newname = tax_vector[4]
      } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
        newname = tax_vector[3]
      } else {
        newname = tax_vector[2]
      }
      newnames = c(newnames,newname)
      top10_otu_table = rbind(top10_otu_table,dd[rownumber,])
    }
    rownames(top10_otu_table)<-newnames
    top10=top10_otu_table
    top10 = top10[,!colnames(top10)=="sum"]
    top10$species<-row.names(top10)
    melt_top10<-melt(top10)
    melt_top10$variable <- factor(melt_top10$variable,levels = cluster_order)
    names(melt_top10)<-c("genus", "ID", "percent")
    p <- plot_ly(data = melt_top10[which(melt_top10$genus == rownames(top10)[10]),], x = ~ID, y = ~percent, type = 'bar', name =rownames(top10)[10]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[9]),], name =rownames(top10)[9]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[8]),], name =rownames(top10)[8]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[7]),], name =rownames(top10)[7]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[6]),], name =rownames(top10)[6]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[5]),], name =rownames(top10)[5]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[4]),], name =rownames(top10)[4]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[3]),], name =rownames(top10)[3]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[2]),], name =rownames(top10)[2]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[1]),], name =rownames(top10)[1]) %>%
      layout(title = plot_title,
             xaxis=list(title="Sample ID"),
             yaxis=list(title="Abundance (percent of rarefied counts)"),
             barmode = 'stack',
             #autosize = F,
             margin = list(l=50,r=50,b=100,t=50))
    p
    tile_df = data.frame("IDs" = sample_names(po), "var_vec" =  split_variable_vector, "y" = rep(1,length(sample_names(po))))
    p3 <- ggplot(tile_df, aes(IDs,y)) +
      geom_tile(aes(fill = var_vec), colour = "grey50") + theme_bw() + labs(fill=test_leg) + scale_fill_manual(values=color_vector)
    p3
    r <- data.frame(ID=sample_names(po), type=factor(get_variable(po,variable_name)), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Shannon")))
    r_ordered = r[order(fit$order),]
    r_ordered$ID = factor(r_ordered$ID, levels = cluster_order)
    rownames(r_ordered) = r_ordered$ID
    p2 <- plot_ly(type = "bar", data = r_ordered, x = ~ID, y = ~Shannon, color = I("#555555")) %>%
      layout(xaxis= list(showticklabels = FALSE))
    return(list(p,p2,p3,fit,r,r_ordered))
  }
}
              
              
set_api_key <- function() {
  Sys.setenv("plotly_username"="thej-ssi")
  Sys.setenv("plotly_api_key"="gFKgrgfaQjKs1GanZdA7")
}

prune_by_variable <- function(po,variable_name,variable_value) {
  return_po = prune_samples(sample_names(po)[which(get_variable(po,variable_name) %in% variable_value)],po)
}

make_PCOA_plot <- function(po,plotname="PCoA_plot") {
  ord <- ordinate(po, method = "PCoA", distance = "bray")
  groups = levels(sample_data(po)$sample_groups)
  if (length(groups) <= 9) {
    col_vec = RColorBrewer::brewer.pal(length(groups),"Set1")
    p = plot_ordination(po,ord, color="sample_groups", title = plotname) + geom_point(size=2, alpha=0.01)+ stat_ellipse(level=0.75) + scale_colour_manual(values = col_vec)
  } else {
    p = plot_ordination(po,ord, color="sample_groups", title = plotname) + geom_point(size=2, alpha=0.01)+ stat_ellipse(level=0.75)
  }
  anosim_test = anosim(t(otu_table(po)),grouping = as.factor(sample_data(po)$sample_groups))
  returnlist = list(p,anosim_test$statistic,anosim_test$signif)
  return(returnlist)	
}

make_PCoA_object <- function(po,variable_name,plot_title="PCoA_plot",color_list=c(),perform_anosim = TRUE,dist_method = "bray",rngseed = 1) {
  set.seed(rngseed)
  if (dist_method == "jaccard" | dist_method == "binary") {
    ord <- ordinate(po, method = "PCoA", distance = "jaccard", binary = TRUE)
  } else {
    ord <- ordinate(po, method = "PCoA", distance = dist_method)
  }
  if (class(get_variable(po,variable_name))=="factor") {
    groups = levels(get_variable(po,variable_name))
  } else {
    groups = levels(factor(get_variable(po,variable_name)))
  }
  if (length(groups) == length(color_list)) {
    col_vec = color_list
    names(col_vec) = groups
    p = plot_ordination(po,ord, color=variable_name, title = plot_title) + geom_point(size=2, alpha=0.01)+ stat_ellipse(level=0.75) + scale_colour_manual(values = col_vec)
  } else if (length(groups) <= 9) {
    print(paste0('Number of colors given (', length(color_list) , ') does not match number of levels in variable (', length(groups),')'))
    col_vec = RColorBrewer::brewer.pal(9,"Set1")[1:length(groups)]
    p = plot_ordination(po,ord, color=as.character(variable_name), title = plot_title) + geom_point(size=2, alpha=0.01)+ stat_ellipse(level=0.75) + scale_colour_manual(values = col_vec)
  } else {
    print(paste0('Number of colors given (', length(color_list) , ') does not match number of levels in variable (', length(groups),')'))
    p = plot_ordination(po,ord, color=as.character(variable_name), title = plot_title) + geom_point(size=2, alpha=0.01)+ stat_ellipse(level=0.75)
  }
  if (perform_anosim) {
    phen_vec = as.character(get_variable(po,variable_name))
    phen_factor = get_variable(po,variable_name)[which(!is.na(phen_vec))]
    if (dist_method == "jaccard" | dist_method == "binary") {
      anosim_test = anosim(t(otu_table(po)[,which(!is.na(phen_vec))]),grouping = phen_factor,permutations = 1000, distance = "jaccard")
    } else {
      anosim_test = anosim(t(otu_table(po)[,which(!is.na(phen_vec))]),grouping = phen_factor,permutations = 1000, distance = dist_method)
    }
    
    returnlist = list("Plot"=p,"Anosim_results"=anosim_test,"Anosim_R"=anosim_test$statistic,"Anosim_p"=anosim_test$signif)
  } else {
    returnlist = list(p)
  }
  return(returnlist)	
}


get_pairwise_distance <- function(dist_mat,ID1, ID2 = NA) {
  if (is.na(ID2)) {
    sub_dist_mat = dist_mat[which(rownames(dist_mat) %in% ID1),which(colnames(dist_mat) %in% ID1)]
    dist_vec = sub_dist_mat[upper.tri(sub_dist_mat)]
  } else {
    dist_vec = as.vector(dist_mat[which(rownames(dist_mat) %in% ID1),which(colnames(dist_mat) %in% ID2)])
  }
  return(dist_vec)
}


extended_anosim <- function(po,variable_name, dist_method = "bray",color_list=c()) {
  if (dist_method == "euclidean") {
    BC_dist = dist(t(otu_table(po)))
  } else {
    BC_dist = vegdist(t(otu_table(po)))
  }
  BC_dist_mat = as.matrix(BC_dist)
  group_factor = get_variable(po,variable_name)
  var_type = as.vector(class(group_factor))
  if (var_type[length(var_type)]=="factor") {
    groups = as.vector(levels(group_factor))
  } else {
    groups = unique(group_factor)
  }
  group_vec = as.vector(group_factor)
  color_vec = setup_color_vector_2(po,variable_name,color_list)[[1]]
  mat = matrix(nrow=0,ncol=3)
  mean_mat = matrix(nrow=0,ncol=4)
  groups_2 = c()
  for (n in 1:length(groups)) {
    for (m in n:length(groups)) {
      ID1 = sample_names(po)[which(group_vec==groups[n])]
      if (n == m) {
        dist_vec = get_pairwise_distance(BC_dist_mat,ID1)
        group_type_vec = rep("Within group",length(dist_vec))
      } else {
        ID2 = sample_names(po)[which(group_vec==groups[m])]
        dist_vec = get_pairwise_distance(BC_dist_mat,ID1,ID2)
        group_type_vec = rep("Between groups",length(dist_vec))
      }
      group_2_name = paste0(groups[n],'-',groups[m])
      print(group_2_name)
      groups_2 = c(groups_2,group_2_name)
      mat_append = cbind(rep(group_2_name,length(dist_vec)),group_type_vec,dist_vec)
      mat = rbind(mat,mat_append)
      mean_mat = rbind(mean_mat,c(groups[n],groups[m],mean(dist_vec),sd(dist_vec)))
    }
  }
  df = as.data.frame(mat)
  colnames(df) = c("Group","Group_type","Distance")
  df$Distance = as.numeric(as.vector(df$Distance))
  p_mat = matrix(nrow=length(groups_2),ncol=length(groups_2))
  for (n in 1:length(groups_2)) {
    for (m in n:length(groups_2)) {
      pval = wilcox.test(df$Distance[which(df$Group==groups_2[n])],df$Distance[which(df$Group==groups_2[m])])$p.value
      p_mat[n,m] = pval
      p_mat[m,n] = pval
    }
  }
  mean_df = as.data.frame(mean_mat)
  colnames(mean_df) = c("Group_1","Group_2","Mean","SD")
  p_df = as.data.frame(p_mat)
  colnames(p_df) = groups_2
  rownames(p_df) = groups_2
  
  
  return(list("Distances"=df, "Means"=mean_df, "p.values"=p_df))
}



check_counts <- function(po) {
  print("Number of sequences found in samples")
  print(sort(colSums(otu_table(po))))
  plot(sort(colSums(otu_table(po))),main = "Number of sequences found in samples")
}



genus_comparison <- function(po,po_split) {
  mat = matrix(nrow = 0, ncol = 2*length(metadata_split_groups)+3)
  for (i in 1:nrow(tax_table(po))) {
    print_vec = c()
    count_list = list()
    for (groupno in 1:length(metadata_split_groups)) {
      group = metadata_split_groups[groupno]
      count_vec = as.vector(otu_table(po_split[[groupno]])[i,])
      mean = mean(count_vec)
      sd = sd(count_vec)
      print_vec = c(print_vec,mean,sd)
      count_list[[groupno]] = count_vec
    }
    count_vec = as.vector(otu_table(po)[i,])
    mean = mean(count_vec)
    sd = sd(count_vec)
    kruskal = kruskal.test(count_list)
    p_value = kruskal$p.value
    print_vec = c(print_vec,mean,sd,p_value)
    mat = rbind(mat,print_vec)
    
  }
  mat = cbind(tax_table(po)[,1:6],mat)
  
  column_names = c()
  for (i in metadata_split_groups) {
    column_names = c(column_names,paste0('Mean - ',i),paste0('SD - ',i))
  }
  
  colnames(mat) = c(colnames(tax_table(rare_prokaryot_genus))[1:6],column_names,'Mean - all','SD - all','p-value')
  df = as.data.frame(mat)
  df$'adjusted p-value' = p.adjust(df$'p-value',method = "bonferroni")
  return(df)
}


taxa_comparison <- function(po,po_split) {
  mat = matrix(nrow = 0, ncol = 2*length(metadata_split_groups)+3)
  for (i in 1:nrow(tax_table(po))) {
    print_vec = c()
    count_list = list()
    for (groupno in 1:length(metadata_split_groups)) {
      group = metadata_split_groups[groupno]
      count_vec = as.vector(otu_table(po_split[[groupno]])[i,])
      mean = mean(count_vec)
      sd = sd(count_vec)
      print_vec = c(print_vec,mean,sd)
      count_list[[groupno]] = count_vec
    }
    count_vec = as.vector(otu_table(po)[i,])
    mean = mean(count_vec)
    sd = sd(count_vec)
    kruskal = kruskal.test(count_list)
    p_value = kruskal$p.value
    print_vec = c(print_vec,mean,sd,p_value)
    mat = rbind(mat,print_vec)
    
  }
  mat = cbind(tax_table(po),mat)
  
  column_names = c()
  for (i in metadata_split_groups) {
    column_names = c(column_names,paste0('Mean - ',i),paste0('SD - ',i))
  }
  
  colnames(mat) = c(colnames(tax_table(po)),column_names,'Mean - all','SD - all','p-value')
  df = as.data.frame(mat)
  df$'adjusted p-value' = p.adjust(df$'p-value',method = "bonferroni")
  return(df)
}


top_taxa_heatmap <- function(po,top_10_taxa,filename) {
  groups = levels(sample_data(po)$sample_groups)
  group_colors = RColorBrewer::brewer.pal(n=length(groups),name="Set1")
  color_names = c("Red","Blue","Green","Purple","Orange","Yellow","Brown","Pink","Grey")
  group_color_vector = as.vector(sample_data(po)$sample_groups)
  
  color_table = matrix(ncol=2,nrow=0)
  for (i in 1:length(groups)) {
    group_color_vector[group_color_vector==groups[i]] = group_colors[i]
    color_table= rbind(color_table,c(groups[i],color_names[i]))
  }
  
  image_filename = paste0(filename,".png")
  legend_filename = paste0(filename,"__color_key.txt")
  write.table(color_table,file = legend_filename,sep="\t",quote = F,col.names = F,row.names = F)
  
  newnames = c()
  taxmat = tax_table(po)
  for (i in 1:length(top_10_taxa)) {
    rownumber = top_10_taxa[i]
    tax_vector = as.vector(taxmat[rownumber,])
    if (!is.na(tax_vector[7])) {
      newname = paste0(tax_vector[6],' ',tax_vector[7])
    } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
      newname = tax_vector[6]
    } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
      newname = tax_vector[5]
    } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
      newname = tax_vector[4]
    } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
      newname = tax_vector[3]
    } else {
      newname = tax_vector[2]
    }
    newnames = c(newnames,newname)
  }
  
  if (length(top_10_taxa) > 10) {
    other_sum = colSums(otu_table(po)[-top_10_taxa,])
    top_10_matrix = rbind(otu_table(po)[top_10_taxa,],other_sum)
    rownames(top_10_matrix) = c(newnames,"Other")
  } else {
    top_10_matrix = otu_table(po)[top_10_taxa,]
    rownames(top_10_matrix) = newnames
  }
  png(filename = image_filename, width = 1000, height = 700, res=72)
  heatmap.2(top_10_matrix,
            trace = "none",
            scale = "row",
            col = colorspace::diverge_hsv(50),
            margins = c(5,15),
            ColSideColors = group_color_vector,
            #main = "Heatmap showing relative abundance of top 10 genera across all samples",
            key.title = "",
            #lwid = 2,
            labCol = F)
  dev.off()
  heatmap.2(top_10_matrix,
            trace = "none",
            scale = "row",
            col = colorspace::diverge_hsv(50),
            margins = c(5,15),
            ColSideColors = group_color_vector,
            #main = "Heatmap showing relative abundance of top 10 genera across all samples",
            key.title = "",
            #lwid = 2,
            labCol = F)
}

make_barplot_plus_object <- function(po,taxa,variable_name,plot_title="",color_vector=c()) {
  if (class(po)=="phyloseq") {
    taxmat = tax_table(po)
    dd<-otu_table(po)
    dd<-apply(dd, 2, function(x) x/sum(x)*100)
    dd<-as.data.frame(dd)
    split_variable_vector = get_variable(po,variable_name)
    dd$sum<-apply(dd, 1, sum)
    dd_sorted<-dd[dd$sum>0,]
    dd_sorted<-dd_sorted[,-ncol(dd_sorted)]
    dd_sorted_d<-vegdist(t(dd_sorted), method="bray")
    fit <- hclust(dd_sorted_d, method="ward.D")
    #plot(fit, cex=0.5) # display dendogram
    cluster_order<-fit$labels[fit$order]
    newnames = c()
    top10_otu_table = matrix(nrow=0,ncol=ncol(dd_sorted))
    for (i in 1:length(taxa)) {
      rownumber = taxa[i]
      tax_vector = as.vector(taxmat[rownumber,])
      if (!is.na(tax_vector[7])) {
        newname = paste0(tax_vector[6],' ',tax_vector[7])
      } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
        newname = tax_vector[6]
      } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
        newname = tax_vector[5]
      } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
        newname = tax_vector[4]
      } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
        newname = tax_vector[3]
      } else {
        newname = tax_vector[2]
      }
      newnames = c(newnames,newname)
      top10_otu_table = rbind(top10_otu_table,dd[rownumber,])
    }
    rownames(top10_otu_table)<-newnames
    top10=top10_otu_table
    top10 = top10[,!colnames(top10)=="sum"]
    top10$species<-row.names(top10)
    melt_top10<-melt(top10)
    melt_top10$variable <- factor(melt_top10$variable,levels = cluster_order)
    names(melt_top10)<-c("genus", "ID", "percent")
    p <- plot_ly(data = melt_top10[which(melt_top10$genus == rownames(top10)[10]),], x = ~ID, y = ~percent, type = 'bar', name =rownames(top10)[10]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[9]),], name =rownames(top10)[9]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[8]),], name =rownames(top10)[8]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[7]),], name =rownames(top10)[7]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[6]),], name =rownames(top10)[6]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[5]),], name =rownames(top10)[5]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[4]),], name =rownames(top10)[4]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[3]),], name =rownames(top10)[3]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[2]),], name =rownames(top10)[2]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[1]),], name =rownames(top10)[1]) %>%
      layout(title = plot_title,
             xaxis=list(title="Sample ID"),
             yaxis=list(title="Abundance (percent of rarefied counts)"),
             barmode = 'stack',
             #autosize = F,
             margin = list(l=50,r=50,b=100,t=50))
    p
    heatmap_data = as.matrix(dd_sorted[1:2,fit$order])
    sorted_variable_vector = as.vector(split_variable_vector)[fit$order]
    variable_levels = levels(split_variable_vector)
    variable_n = length(variable_levels)
    if (missing(color_vector) | !length(color_vector)==variable_n) {
      if (variable_n < 10) {
        Rcol_vec = RColorBrewer::brewer.pal(variable_n,"Set1")
      } else if (variable_n < 13) {
        Rcol_vec = RColorBrewer::brewer.pal(variable_n,"Set3")
      } else {
        Rcol_vec = grDevices::rainbow(variable_n)
      }
    } else {
      Rcol_vec = color_vector
    }
    
    col_vec = sorted_variable_vector
    for (i in 1:length(Rcol_vec)) {
      col_vec = replace(col_vec,col_vec==variable_levels[i],Rcol_vec[i])
    }
    heatmap.2(heatmap_data,
              Rowv = FALSE,
              Colv = FALSE,
              trace = "none",
              scale = "row",
              dendrogram = "none",
              col = colorspace::diverge_hsv(50),
              margins = c(5,15),
              ColSideColors = col_vec,
              #main = "Heatmap showing relative abundance of top 10 genera across all samples",
              key.title = "")
    r <- data.frame(ID=sample_names(po), type=factor(get_variable(po,variable_name)), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Shannon")))
    r_ordered = r[order(fit$order),]
    r_ordered$ID = factor(r_ordered$ID, levels = cluster_order)
    rownames(r_ordered) = r_ordered$ID
    p2 <- plot_ly(type = "bar", data = r_ordered, x = ~ID, y = ~Shannon, color = I("#555555")) %>%
      layout(xaxis= list(showticklabels = FALSE))
    return(list(p,p2,fit,r,r_ordered))
  }
}
              
              
make_ordered_barplot <- function(po,taxa,plot_variable,sort_variable,plot_title="",col_vec=c()) {
  if (class(po)=="phyloseq") {
    color_vector = setup_color_vector(po,plot_variable,col_vec)
    taxmat = tax_table(po)
    dd<-otu_table(po)
    dd<-apply(dd, 2, function(x) x/sum(x)*100)
    dd<-as.data.frame(dd)
    split_variable_vector = get_variable(po,plot_variable)
    dd$sum<-apply(dd, 1, sum)
    dd_sorted<-dd[dd$sum>0,]
    dd_sorted<-dd_sorted[,-ncol(dd_sorted)]
    #plot(fit, cex=0.5) # display dendogram
    newnames = c()
    top10_otu_table = matrix(nrow=0,ncol=ncol(dd_sorted))
    for (i in 1:length(taxa)) {
      rownumber = taxa[i]
      tax_vector = as.vector(taxmat[rownumber,])
      if (!is.na(tax_vector[7])) {
        newname = paste0(tax_vector[6],' ',tax_vector[7])
      } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
        newname = tax_vector[6]
      } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
        newname = tax_vector[5]
      } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
        newname = tax_vector[4]
      } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
        newname = tax_vector[3]
      } else {
        newname = tax_vector[2]
      }
      newnames = c(newnames,newname)
      top10_otu_table = rbind(top10_otu_table,dd[rownumber,])
    }
    rownames(top10_otu_table)<-newnames
    top10=top10_otu_table
    top10 = top10[,!colnames(top10)=="sum"]
    top10$species<-row.names(top10)
    melt_top10<-melt(top10)
    melt_top10$variable <- factor(melt_top10$variable,levels = sample_names(po)[order(get_variable(po,sort_variable))])
    names(melt_top10)<-c("genus", "ID", "percent")
    p <- plot_ly(data = melt_top10[which(melt_top10$genus == rownames(top10)[10]),], x = ~ID, y = ~percent, type = 'bar', name =rownames(top10)[10]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[9]),], name =rownames(top10)[9]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[8]),], name =rownames(top10)[8]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[7]),], name =rownames(top10)[7]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[6]),], name =rownames(top10)[6]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[5]),], name =rownames(top10)[5]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[4]),], name =rownames(top10)[4]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[3]),], name =rownames(top10)[3]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[2]),], name =rownames(top10)[2]) %>%
      add_trace(data = melt_top10[which(melt_top10$genus == rownames(top10)[1]),], name =rownames(top10)[1]) %>%
      layout(title = plot_title,
             xaxis=list(title="Sample ID"),
             yaxis=list(title="Abundance (percent of rarefied counts)"),
             barmode = 'stack',
             #autosize = F,
             margin = list(l=50,r=50,b=100,t=50))
    p
    tile_df = data.frame("Type"=as.character(as.vector(get_variable(po,plot_variable)))[order(get_variable(po,sort_variable))],"y" = rep(1,length(sample_names(po))),"ID"=sample_names(po)[order(get_variable(po,sort_variable))])
    p_tiles = ggplot(data=tile_df)+geom_tile(aes(x=ID, y=y, fill=Type))+theme_bw()+scale_fill_manual(values=color_vector)+theme(axis.text.x = element_text(angle = -90, hjust = 1))
    ggplot(data=tile_df)+geom_tile(aes(x=ID, y=y, fill=Type))+theme_bw()+scale_fill_manual(values=col_vec)+theme(axis.text.x = element_text(angle = -90, hjust = 1))
    r <- data.frame(ID=sample_names(po), type=factor(get_variable(po,plot_variable)), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Shannon")))
    r_ordered = r[order(get_variable(po,sort_variable)),]
    r_ordered$ID = factor(r_ordered$ID, levels = sample_names(po)[order(get_variable(po,sort_variable))])
    rownames(r_ordered) = r_ordered$ID
    p2 <- plot_ly(type = "bar", data = r_ordered, x = ~ID, y = ~Shannon, color = I("#555555")) %>%
      layout(xaxis= list(showticklabels = FALSE))
    return(list(p,p_tiles))
  }
}

make_heatmap_object <- function(po,top_10_taxa,variable_name,plot_title="Heatmap",color_list=c()) {
  group_color_vector = as.vector(get_variable(po,variable_name))
  groups = levels(factor(group_color_vector))
  group_count = length(groups)
  if (length(color_list) == group_count) {
    group_colors = color_list
  } else if (group_count<=9) {
    group_colors = RColorBrewer::brewer.pal(n=length(groups),name="Set1")
  } else {
    group_colors = grDevices::rainbow(group_count)
  }
  color_table = matrix(ncol=2,nrow=0)
  for (i in 1:length(groups)) {
    group_color_vector[group_color_vector==groups[i]] = group_colors[i]
    color_table= rbind(color_table,c(groups[i],group_colors[i]))
  }
  newnames = c()
  taxmat = tax_table(po)
  for (i in 1:length(top_10_taxa)) {
    rownumber = top_10_taxa[i]
    tax_vector = as.vector(taxmat[rownumber,])
    if (!is.na(tax_vector[7])) {
      newname = paste0(tax_vector[6],' ',tax_vector[7])
    } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
      newname = tax_vector[6]
    } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
      newname = tax_vector[5]
    } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
      newname = tax_vector[4]
    } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
      newname = tax_vector[3]
    } else {
      newname = tax_vector[2]
    }
    newnames = c(newnames,newname)
  }
  
  if (length(top_10_taxa) > 10) {
    other_sum = colSums(otu_table(po)[-top_10_taxa,])
    top_10_matrix = rbind(otu_table(po)[top_10_taxa,],other_sum)
    rownames(top_10_matrix) = c(newnames,"Other")
  } else {
    top_10_matrix = otu_table(po)[top_10_taxa,]
    rownames(top_10_matrix) = newnames
  }
  print(group_color_vector)
  return_object = heatmap.2(top_10_matrix,
                            trace = "none",
                            scale = "row",
                            col = colorspace::diverge_hsv(50),
                            margins = c(5,15),
                            ColSideColors = group_color_vector,
                            #main = "Heatmap showing relative abundance of top 10 genera across all samples",
                            key.title = "",
                            #lwid = 2,
                            labCol = F)
  print(color_table)
  return(return_object)
}



make_violin_object <- function(po,variable_name,taxa,plot_title="Distribution of relative abundance",color_list=c(),level_list=c()) {
  color_vector = setup_color_vector(po,variable_name,color_list)
  otu = otu_table(po)[taxa,]
  tax_vector = get_taxa_names(po,taxa)
  variable_vector = as.vector(get_variable(po,variable_name))
  data = rbind(otu,variable_vector)
  data = as.data.frame(t(data))
  colnames(data) = c(tax_vector,"Group")
  data2 = melt.data.frame(data,id.vars = "Group")
  data2$value = as.numeric(as.vector(data2$value))
  if (!missing(level_list) & length(level_list)==length(levels(factor(data2$Group)))) {
    p <- ggplot(data2, aes(x=variable, y=value, fill = factor(Group,levels = level_list))) + 
      geom_violin() + 
      coord_flip() +
      scale_y_log10() +
      scale_fill_manual(values = color_vector,name = variable_name) +
      labs(y = "Rarefied sequence counts", x = "Taxa", title = plot_title)
    p
  } else {
    p <- ggplot(data2, aes(x=variable, y=value, fill = Group)) + 
      geom_violin() + 
      coord_flip() +
      scale_y_log10() +
      scale_fill_manual(values = color_vector,name = variable_name) +
      labs(y = "Rarefied sequence counts", x = "Taxa", title = plot_title)
    p
  }
  
  return(list(p,data,data2))
}

make_taxa_comparison_object <- function(po,variable_name,p_adjust_method="bonferroni") {
  d = otu_table(po)
  tax = tax_table(po)
  variable_vector = as.vector(get_variable(po,variable_name))
  groups = unique(variable_vector)
  variable_count = length(groups)
  p_mat = tax
  if (variable_count < 2) {
    print("Less than two types found in designated variable")
  } else if (variable_count == 2) {
    mean_headers = c('mean, all')
    group_means = cbind(apply(d,1, function(e) mean(as.numeric(e))))
    for (group in groups) {
      group_mean = apply(d,1, function(e) mean(as.numeric(e[which(variable_vector == group)])))
      group_means = cbind(group_means,group_mean)
      mean_headers = c(mean_headers,paste0('mean, ',group))
    }
    MWU_tests = apply(d,1, function(e) wilcox.test(as.numeric(e[which(variable_vector == groups[1])]),as.numeric(e[which(variable_vector == groups[2])]))$p.value)
    MWU_corrected = p.adjust(MWU_tests,method = p_adjust_method)
    mwu_headers = c('p','p.adjusted')
    p_mat = cbind(p_mat,as.numeric(MWU_tests),as.numeric(MWU_corrected))
  } else {
    mean_headers = c('mean, all')
    group_means = cbind(apply(d,1, function(e) mean(as.numeric(e))))
    for (group in groups) {
      group_mean = apply(d,1, function(e) mean(as.numeric(e[which(variable_vector == group)])))
      group_means = cbind(group_means,group_mean)
      mean_headers = c(mean_headers,paste0('mean, ',group))
    }
    mwu_headers = c()
    for (n1 in 1:(length(groups)-1)) {
      group1 = groups[n1]
      for (n2 in (n1+1):length(groups)) {
        group2 = groups[n2]
        MWU_tests = apply(d,1, function(e) wilcox.test(as.numeric(e[which(variable_vector == group1)]),as.numeric(e[which(variable_vector == group2)]))$p.value)
        MWU_corrected = p.adjust(MWU_tests,method = p_adjust_method)
        mwu_header = paste0('p, ',group1,' v ',group2)
        mwu_corrected_header = paste0('p.adjusted_',group1,'.v.',group2)
        mwu_headers = c(mwu_headers,mwu_header,mwu_corrected_header)
        p_mat = cbind(p_mat,as.numeric(MWU_tests),as.numeric(MWU_corrected))
      }
    }
  }
  return_df = as.data.frame(cbind(rep(1:nrow(p_mat)),p_mat,group_means))
  #apply(return_df[,9:ncol(return_df)],2, function(e) as.numeric(e))
  colnames(return_df) = c("Rownumber",colnames(tax),mwu_headers,mean_headers)
  rownames(return_df) = rownames(tax)
  return(return_df)
}
                          

make_paired_taxa_comparison <- function(po,variable_name,ID_name,p_adjust_method="bonferroni") {
  d = otu_table(po)
  tax = tax_table(po)
  variable_factor = get_variable(po,variable_name)
  variable_vector = as.vector(variable_factor)
  if (class(variable_factor)=="factor") {
    groups = levels(variable_factor)
  } else {
    groups = unique(variable_vector)
  }
  ID_vector = as.vector(get_variable(po,ID_name))
  ID_1_rows = which(variable_vector == groups[1])
  ID_2_rows = which(variable_vector == groups[2])
  ID_vector_1 = ID_vector[ID_1_rows]
  ID_vector_2 = ID_vector[ID_2_rows]
  ID_1_rows = ID_1_rows[which(ID_vector_1 %in% ID_vector_2)]
  ID_2_rows = ID_2_rows[which(ID_vector_2 %in% ID_vector_1)]
  include_index = c(ID_1_rows[order(ID_vector_1)],ID_2_rows[order(ID_vector_2)])
  print(include_index)
  d = d[,include_index]
  print(variable_vector)
  variable_vector = variable_vector[include_index]
  print(variable_vector)
  variable_count = length(groups)
  p_mat = tax
  if (variable_count < 2) {
    print("Less than two types found in designated variable")
  } else if (variable_count == 2) {
    mean_headers = c('mean, all')
    
    group_means = cbind(apply(d,1, function(e) mean(as.numeric(e))))
    for (group in groups) {
      group_mean = apply(d,1, function(e) mean(as.numeric(e[which(variable_vector == group)])))
      group_means = cbind(group_means,group_mean)
      mean_headers = c(mean_headers,paste0('mean, ',group))
    }
    MWU_tests = apply(d,1, function(e) wilcox.test(as.numeric(e[which(variable_vector == groups[1])]),as.numeric(e[which(variable_vector == groups[2])]), paired = TRUE)$p.value)
    MWU_corrected = p.adjust(MWU_tests,method = p_adjust_method)
    mwu_headers = c('p','p.adjusted')
    p_mat = cbind(p_mat,as.numeric(MWU_tests),as.numeric(MWU_corrected))
   }
  #  else {
  #   mean_headers = c('mean, all')
  #   group_means = cbind(apply(d,1, function(e) mean(as.numeric(e))))
  #   for (group in groups) {
  #     group_mean = apply(d,1, function(e) mean(as.numeric(e[which(variable_vector == group)])))
  #     group_means = cbind(group_means,group_mean)
  #     mean_headers = c(mean_headers,paste0('mean, ',group))
  #   }
  #   mwu_headers = c()
  #   for (n1 in 1:(length(groups)-1)) {
  #     group1 = groups[n1]
  #     for (n2 in (n1+1):length(groups)) {
  #       group2 = groups[n2]
  #       MWU_tests = apply(d,1, function(e) wilcox.test(as.numeric(e[which(variable_vector == group1)]),as.numeric(e[which(variable_vector == group2)]))$p.value)
  #       MWU_corrected = p.adjust(MWU_tests,method = p_adjust_method)
  #       mwu_header = paste0('p, ',group1,' v ',group2)
  #       mwu_corrected_header = paste0('p.adjusted_',group1,'.v.',group2)
  #       mwu_headers = c(mwu_headers,mwu_header,mwu_corrected_header)
  #       p_mat = cbind(p_mat,as.numeric(MWU_tests),as.numeric(MWU_corrected))
  #     }
  #   }
  # }
  return_df = as.data.frame(cbind(rep(1:nrow(p_mat)),p_mat,group_means))
  #apply(return_df[,9:ncol(return_df)],2, function(e) as.numeric(e))
  colnames(return_df) = c("Rownumber",colnames(tax),mwu_headers,mean_headers)
  rownames(return_df) = rownames(tax)
  if ('p' %in% colnames(return_df)) {
    return_df$p = as.numeric(as.vector(return_df$p))
    return_df$p.adjusted = as.numeric(as.vector(return_df$p.adjusted))
  }
  return(return_df)
}                         
                          
make_presence_absence_comparison_object <- function(po,variable_name,p_adjust_method="bonferroni",presence_threshold_percent = 0) {
  d = otu_table(po)
  tax = tax_table(po)
  variable_vector = as.vector(get_variable(po,variable_name))
  groups = unique(variable_vector)
  variable_count = length(groups)
  p_mat = tax
  if (variable_count < 2) {
    print("Less than two types found in designated variable")
  } else {
    d2 = d
    d2[d<=presence_threshold_percent] <- 0
    d2[d>presence_threshold_percent] <- 1
    d_g1 = d2[,which(variable_vector == groups[1])]
    d_g2 = d2[,which(variable_vector == groups[2])]
    sample_count = ncol(d)
    sample_count_g1 = ncol(d_g1)
    sample_count_g2 = ncol(d_g2)
    pval <- apply(d2, 1, function(e) fisher.test(table(factor(e, levels=c(0,1)), variable_vector))$p.value)
    g1_presence_vector <- apply(d_g1, 1, sum)
    g1_absence_vector = sample_count-g1_presence_vector
    g2_presence_vector <- apply(d_g2, 1, sum)
    g2_absence_vector = sample_count-g2_presence_vector
    p.corrected = p.adjust(pval,method = p_adjust_method)
    return_df = as.data.frame(cbind(rep(1:nrow(p_mat)),tax,g1_presence_vector,g1_absence_vector,g2_presence_vector,g2_absence_vector,pval,p.corrected))
    colnames(return_df) = c("Rownumber",colnames(tax),paste0(groups[1],'_present'),paste0(groups[1],'_absent'),paste0(groups[2],'_present'),paste0(groups[2],'_absent'),"p","p.corrected")
    rownames(return_df) = rownames(tax)
    return(return_df)
  }
}


make_OTU_boxplot_object_2 <- function(po,OTU,variable_name,plot_name="Distribution of relative abundance",color_list=c()) {
  variable_vector = as.vector(get_variable(po,variable_name))
  groups = unique(variable_vector)
  d = as.vector(otu_table(po)[OTU,])
  tax_vector = as.vector(tax_table(po)[OTU,])
  if (!is.na(tax_vector[7])) {
    newname = paste0(tax_vector[6],' ',tax_vector[7])
  } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
    newname = tax_vector[6]
  } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
    newname = tax_vector[5]
  } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
    newname = tax_vector[4]
  } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
    newname = tax_vector[3]
  } else {
    newname = tax_vector[2]
  }
  col_vec = setup_color_vector(po,variable_name,color_list)
  p <- plot_ly(y = d, color = variable_vector, type = "box", boxpoints = "all", pointpos = -1.5, colors = col_vec) %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = 90),
           yaxis=list(title=paste0(newname, ' rarefied sequence counts')),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
}
                          

make_OTU_boxplot_object <- function(po,OTU,variable_name,plot_name="",color_list=c()) {
  variable_factor = get_variable(po,variable_name)
  groups = levels(variable_factor)
  d = as.vector(otu_table(po)[OTU,])
  tax_vector = as.vector(tax_table(po)[OTU,])
  if (!is.na(tax_vector[7])) {
    newname = paste0(tax_vector[6],' ',tax_vector[7])
  } else if (!is.na(tax_vector[6]) & !tax_vector[6]=="unclassified") {
    newname = tax_vector[6]
  } else if (!is.na(tax_vector[5]) & !tax_vector[5]=="unclassified") {
    newname = tax_vector[5]
  } else if (!is.na(tax_vector[4]) & !tax_vector[4]=="unclassified") {
    newname = tax_vector[4]
  } else if (!is.na(tax_vector[3]) & !tax_vector[3]=="unclassified") {
    newname = tax_vector[3]
  } else {
    newname = tax_vector[2]
  }
  col_vec = setup_color_vector(po,variable_name,color_list)
  p <- plot_ly(y = d, color = variable_factor, type = "box", boxpoints = "all", pointpos = -1.5, colors = col_vec) %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = 90),
           yaxis=list(title=paste0(newname, ' rarefied sequence counts')),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  p_mat = matrix(ncol = length(groups), nrow = length(groups))
  for (n1 in 1:(length(groups)-1)) {
    print(n1)
    for (n2 in (n1+1):length(groups)) {
      group1 = groups[n1]
      group2 = groups[n2]
      print(paste0('group1 ',group1))
      print(paste0('group2 ',group2))
      vec1 = d[which(variable_factor==group1)]
      vec2 = d[which(variable_factor==group2)]
      wilcox_test = wilcox.test(vec1,vec2)
      p_mat[n1,n2] <- wilcox_test$p.value
      p_mat[n2,n1] <- wilcox_test$p.value
    }
    
  }
  p_df = as.data.frame(p_mat)
  rownames(p_df) = groups
  colnames(p_df) = groups
  return_list = list('plot'=p,'p.values'=p_df)
}




setup_color_vector <- function(po,variable_name,color_list) {
  group_color_vector = get_variable(po,variable_name)
  groups = levels(group_color_vector)
  group_count = length(groups)
  if (length(color_list) == group_count) {
    group_colors = color_list
  } else if (group_count<=9) {
    group_colors = RColorBrewer::brewer.pal(9,name="Set1")[1:group_count]
  } else if (group_count<=12) {
    group_colors = RColorBrewer::brewer.pal(12,name="Set3")[1:group_count]
  } else {
    group_colors = grDevices::rainbow(group_count)
  }
  return(group_colors)
}
                          
setup_color_vector_2 <- function(po,variable_name,color_list) {
  group_color_vector = as.vector(get_variable(po,variable_name))
  groups = levels(factor(group_color_vector))
  group_count = length(groups)
  if (length(color_list) == group_count) {
    group_colors = color_list
  } else if (group_count<=9) {
    group_colors = RColorBrewer::brewer.pal(9,name="Set1")[1:group_count]
  } else if (group_count<=12) {
    group_colors = RColorBrewer::brewer.pal(12,name="Set3")[1:group_count]
  } else {
    group_colors = grDevices::rainbow(group_count)
  }
  color_table = matrix(ncol=2,nrow=0)
  for (i in 1:length(groups)) {
    group_color_vector[group_color_vector==groups[i]] = group_colors[i]
    color_table= rbind(color_table,c(groups[i],group_colors[i]))
  }
  return(list(group_colors,color_table))
}

                  
                  
test_color_tile <- function(color_vec) {
  values = rep(1,length(color_vec))
  col_vec = factor(color_vec,levels=c(as.character(color_vec)))
  p <- plot_ly(type = "bar", x = col_vec, y = ~values, name = col_vec, color = col_vec, colors = color_vec)
  return(p)
}


make_legend_color <- function(po,variable_name,color_list=c()) {
  if (class(get_variable(po,variable_name))=="factor") {
    groups = levels(get_variable(po,variable_name))
  } else {
    groups = levels(factor(get_variable(po,variable_name)))
  }
  values = rep(1,length(groups))
  variable_n = length(groups)
  Rcol_vec = setup_color_vector(po,variable_name,color_list = color_list)
  xgroups = groups
  p <- plot_ly(type = "bar", x = xgroups, y = values, name = xgroups, color = xgroups, colors = Rcol_vec)
  return(p)
}


plot_colors <- function(color_vec) {
  df = data.frame('ID'=names(color_vec),'color'=color_vec,'count'=1)
  df$color = as.vector(df$color)
  p <- plot_ly(data = df, type='bar',x=~ID,y=~count,color=~ID,colors=~color)
  p
}
               
                  

get_alphadiversity_change_observed <- function(po,pt_var,time_var,time_values,change_var) {
  m = sample_data(po)
  d = otu_table(po)
  otu_change_mat = matrix(nrow = nrow(d),ncol=0)
  m_mat = matrix(nrow=0,ncol=ncol(m))
  patient_vec = get_variable(po,pt_var)
  time_vec = get_variable(po,time_var)
  patients = levels(patient_vec)
  ID_vec = sample_names(po)
  IDs = c()
  change_vec = c()
  change_var_vec = get_variable(po,change_var)
  AD_change_vec = c()
  r <- data.frame(ID=sample_names(po), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon")))
  for (pt in patients) {
    row1 = which(patient_vec == pt & time_vec == time_values[1])
    row2 = which(patient_vec == pt & time_vec == time_values[2])
    if (length(row1) > 0 & length(row2) > 0) {
      AD_change = r$Observed[row2]-r$Observed[row1]
      print(AD_change)
      AD_change_vec = c(AD_change_vec,AD_change)
      IDs = c(IDs,pt)
      change_1 = change_var_vec[row1]
      change_2 = change_var_vec[row2]
      m_mat = rbind(m_mat,m[row1,])
      if (class(change_1)=="numeric" & class(change_2)=="numeric") {
        change = change_2-change_1
      } else {
        change = paste0(change_1,' - ',change_2)
      }
      change_vec = c(change_vec,change)
    }
  }
  m_df = as.data.frame(m_mat)
  rownames(m_df) = IDs
  return(list(AD_change_vec,m_df,change_vec))
}

get_alphadiversity_change_shannon <- function(po,pt_var,time_var,time_values,change_var) {
  m = sample_data(po)
  d = otu_table(po)
  otu_change_mat = matrix(nrow = nrow(d),ncol=0)
  m_mat = matrix(nrow=0,ncol=ncol(m))
  patient_vec = get_variable(po,pt_var)
  time_vec = get_variable(po,time_var)
  patients = levels(patient_vec)
  ID_vec = sample_names(po)
  IDs = c()
  change_vec = c()
  change_var_vec = get_variable(po,change_var)
  AD_change_vec = c()
  r <- data.frame(ID=sample_names(po), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon")))
  for (pt in patients) {
    row1 = which(patient_vec == pt & time_vec == time_values[1])
    row2 = which(patient_vec == pt & time_vec == time_values[2])
    if (length(row1) > 0 & length(row2) > 0) {
      AD_change = r$Shannon[row2]-r$Shannon[row1]
      print(AD_change)
      AD_change_vec = c(AD_change_vec,AD_change)
      IDs = c(IDs,pt)
      change_1 = change_var_vec[row1]
      change_2 = change_var_vec[row2]
      m_mat = rbind(m_mat,m[row1,])
      if (class(change_1)=="numeric" & class(change_2)=="numeric") {
        change = change_2-change_1
      } else {
        change = paste0(change_1,' - ',change_2)
      }
      change_vec = c(change_vec,change)
    }
  }
  m_df = as.data.frame(m_mat)
  rownames(m_df) = IDs
  return(list(AD_change_vec,m_df,change_vec))
}
                  
                  

run_cross_sectional_analysis <- function(po, variable_name, color_list) {
  sample_data(po)$Group = as.character(as.vector(get_variable(po,variable_name)))
  color_vector = setup_color_vector(po,"Group",color_list)
  
  Alphadiv_plot = make_alphadiversity_object(po,variable_name = "Group",plot_title = paste0("Alpha diversity of ",variable_name),color_vector)
  Alphadiv_plot
  
  PCoA_plot = make_PCoA_object(po,variable_name = "Group",plot_title = paste0("Principal coordinate analysis of ",variable_name),color_vector)
  PCoA_plot
  
  top10_taxa = get_top_n_taxa(po,10)
  
  bar_plot = make_abundance_barplot(po,taxa = top10_taxa, "Top 10 most abundant species")
  bar_plot
  make_barplot_plus_object(po,top10_taxa,"Group","Top 10 most abundant species",color_vector)
  
  Heatmap = make_heatmap_object(po,get_top_n_taxa(po,30),"Group",paste0("Heatmap showing over and underrepresentation of top 30 species, ",variable_name),color_vector)
  Heatmap
  
  taxa_comparison_df = make_taxa_comparison_object(po,"Group","bonferroni")
  
  return(list(Alphadiv_plot,PCoA_plot,bar_plot,Heatmap,taxa_comparison_df))
}
                          

get_alphadiversity_change <- function(po,pt_var,time_var,time_values,change_var) {
  m = sample_data(po)
  d = otu_table(po)
  otu_change_mat = matrix(nrow = nrow(d),ncol=0)
  m_mat = matrix(nrow=0,ncol=ncol(m))
  patient_vec = get_variable(po,pt_var)
  time_vec = get_variable(po,time_var)
  if (class(patient_vec) == "factor") {
    patients = levels(patient_vec)
  } else {
    patients = unique(patient_vec)
  }
  ID_vec = sample_names(po)
  IDs = c()
  change_vec = c()
  change_var_vec = get_variable(po,change_var)
  AD_change_vec = c()
  shannon_change_vec = c()
  r <- data.frame(ID=sample_names(po), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon")))
  for (pt in patients) {
    row1 = which(patient_vec == pt & time_vec == time_values[1])
    row2 = which(patient_vec == pt & time_vec == time_values[2])
    if (length(row1) > 0 & length(row2) > 0) {
      AD_change = r$Observed[row2]-r$Observed[row1]
      shannon_change = r$Shannon[row2]-r$Shannon[row1]
      change_1 = change_var_vec[row1]
      change_2 = change_var_vec[row2]
      if (!is.na(change_1) & !is.na(change_2)) {
        IDs = c(IDs,pt)
        m_mat = rbind(m_mat,m[row1,])
        AD_change_vec = c(AD_change_vec,AD_change)
        shannon_change_vec = c(shannon_change_vec,shannon_change)
        if (class(change_1) %in% c("integer","numeric") & class(change_2) %in% c("integer","numeric")) {
          change = change_2-change_1
        } else {
          change = paste0(change_1,' - ',change_2)
        }
        change_vec = c(change_vec,change)
      }
    }
  }
  m_df = as.data.frame(m_mat)
  rownames(m_df) = IDs
  m_df$delta_richness = AD_change_vec
  m_df$delta_Shannon = shannon_change_vec
  m_df$variable_change = change_vec
  return(m_df)
}

get_alphadiversity_change_consecutive <- function(po,pt_var,time_var,time_values,change_var) {
  m_list = list()
  for (i in 1:(length(time_values)-1)) {
    t1 = time_values[i]
    t2 = time_values[i+1]
    sub_m = get_alphadiversity_change(po,pt_var,time_var,c(t1,t2),change_var)
    if (nrow(sub_m) > 0) {
      sub_m$time_change = paste0(t1,' - ',t2)
      m_list[[i]] = sub_m
    }
  }
  m_combined = m_list[[1]]
  for (i in 2:length(m_list)) {
    m_combined = rbind(m_combined,m_list[[i]])
  }
  return(m_combined)
}


get_alphadiversity_change_from_zero <- function(po,pt_var,time_var,time_values,change_var) {
  m_list = list()
  t1 = time_values[1]
  for (i in 1:(length(time_values)-1)) {
    t2 = time_values[i+1]
    sub_m = get_alphadiversity_change(po,pt_var,time_var,c(t1,t2),change_var)
    if (nrow(sub_m) > 0) {
      sub_m$time_change = paste0(t2)
      m_list[[i]] = sub_m
    }
  }
  m_combined = m_list[[1]]
  for (i in 2:length(m_list)) {
    m_combined = rbind(m_combined,m_list[[i]])
  }
  return(m_combined)
}




run_all <- function(rare_prokaryot,rare_eukaryot,rare_fungi,summary_list) {
  ### Setup directory for plots and tables ###
  output_dir = setup_outdir()
  
  ### Create genus level phyloseq objects
  rare_prokaryot_genus<-tax_glom(rare_prokaryot, taxrank = "Genus")
  rare_eukaryot_genus<-tax_glom(rare_eukaryot, taxrank = "Genus")
  rare_fungi_genus<-tax_glom(rare_fungi, taxrank = "Genus")
  
  ### Create phyloseq objects for subsets of data ###
  rare_split_prokaryot_genus = split_by_metadata(rare_prokaryot_genus,metadata_split_variable,metadata_split_groups)
  rare_split_eukaryot_genus = split_by_metadata(rare_eukaryot_genus,metadata_split_variable,metadata_split_groups)
  rare_split_fungi_genus = split_by_metadata(rare_fungi_genus,metadata_split_variable,metadata_split_groups)
  
  
  ### Create alphadiversity plots ###
  alphadiversity_prokaryot = alphadiversity_plot(rare_prokaryot,'Prokaryotic alphadiversity')
  filename = paste0(output_dir,"/Fig_1-",1,"_observed_diversity_prokaryot",".png")
  export(alphadiversity_prokaryot[[1]],file=filename)
  filename = paste0(output_dir,"/Fig_1-",2,"_Shannon_diversity_prokaryot",".png")
  export(alphadiversity_prokaryot[[2]],file=filename)
  
  alphadiversity_eukaryot = alphadiversity_plot(rare_eukaryot,'Eukaryotic alphadiversity')
  filename = paste0(output_dir,"/Fig_1-",3,"_observed_diversity_eukaryot",".png")
  export(alphadiversity_eukaryot[[1]],file=filename)
  filename = paste0(output_dir,"/Fig_1-",4,"_Shannon_diversity_eukaryot",".png")
  export(alphadiversity_eukaryot[[2]],file=filename)
  
  alphadiversity_fungi = alphadiversity_plot(rare_fungi,'Fungal alphadiversity')
  filename = paste0(output_dir,"/Fig_1-",5,"_observed_diversity_fungi",".png")
  export(alphadiversity_fungi[[1]],file=filename)
  filename = paste0(output_dir,"/Fig_1-",6,"_Shannon_diversity_fungi",".png")
  export(alphadiversity_fungi[[2]],file=filename)
  
  
  ### PCoA plots
  PCoA_prokaryot<-make_PCOA_plot(rare_prokaryot,'PCoA based on Bray-curtis dissimilarity of prokaryot taxa')
  ggsave(filename = paste0(output_dir,"/Fig_2-1_PCoA_prokaryot.png"),plot = PCoA_prokaryot[[1]])
  
  PCoA_eukaryot<-make_PCOA_plot(rare_eukaryot,'PCoA based on Bray-curtis dissimilarity of eukaryot taxa')
  ggsave(filename = paste0(output_dir,"/Fig_2-2_PCoA_eukaryot.png"),plot = PCoA_eukaryot[[1]])
  
  PCoA_fungi<-make_PCOA_plot(rare_fungi,'PCoA based on Bray-curtis dissimilarity of fungal taxa')
  ggsave(filename = paste0(output_dir,"/Fig_2-3_PCoA_fungi.png"),plot = PCoA_fungi[[1]])
  
  
  ### get top 10 genus ###
  top10_prokaryot_genus = get_top10_taxa(rare_prokaryot_genus)
  top10_eukaryot_genus = get_top10_taxa(rare_eukaryot_genus)
  top10_fungi_genus = get_top10_taxa(rare_fungi_genus)  
  
  ### Abundance barplots for all samples ###
  p = make_abundance_barplot(rare_prokaryot_genus,top10_prokaryot_genus,'Relative abundance of top 10 prokaryotic genera in all samples')
  print(p)
  filename = paste0(output_dir,"/Fig_3-1_abundance_prokaryot_all.png")
  export(p,file=filename)
  
  p = make_abundance_barplot(rare_eukaryot_genus,top10_eukaryot_genus,'Relative abundance of top 10 eukaryotic genera in all samples')
  print(p)
  filename = paste0(output_dir,"/Fig_3-2_abundance_eukaryot_all.png")
  export(p,file=filename)
  
  p = make_abundance_barplot(rare_fungi_genus,top10_fungi_genus,'Relative abundance of top 10 fungal genera in all samples')
  print(p)
  filename = paste0(output_dir,"/Fig_3-3_abundance_fungi_all.png")
  export(p,file=filename)
  
  
  
  ### Create Abundance barplot for each subset of samples, prokaryot ###
  for (i in 1:length(metadata_split_groups)) {
    plot_name = paste0('Relative abundance of top 10 prokaryotic genera in ',metadata_split_groups[i],' samples')
    p=make_abundance_barplot(rare_split_prokaryot_genus[[i]],top10_prokaryot_genus,plot_name)
    print(p)
    filename = paste0(output_dir,"/Fig_4-",i,"_abundance_prokaryot_",metadata_split_groups[i],".png")
    export(p,file=filename)
  }
  
  ### Create Abundance barplot for each subset of samples, eukaryota ###
  for (i in 1:length(metadata_split_groups)) {
    plot_name = paste0('Relative abundance of top 10 eukaryotic genera in ',metadata_split_groups[i],' samples')
    p=make_abundance_barplot(rare_split_eukaryot_genus[[i]],top10_eukaryot_genus,plot_name)
    print(p)
    filename = paste0(output_dir,"/Fig_5-",i,"_abundance_eukaryot_",metadata_split_groups[i],".png")
    export(p,file=filename)
  }
  
  ### Create Abundance barplot for each subset of samples, fungi ###
  for (i in 1:length(metadata_split_groups)) {
    plot_name = paste0('Relative abundance of top 10 fungal genera in ',metadata_split_groups[i],' samples')
    p=make_abundance_barplot(rare_split_fungi_genus[[i]],top10_fungi_genus,plot_name)
    print(p)
    filename = paste0(output_dir,"/Fig_6-",i,"_abundance_fungi_",metadata_split_groups[i],".png")
    export(p,file=filename)
  }
  
  
  ### Create Abundance barplot for each subset of samples with individual top 10 list, prokaryot###
  for (i in 1:length(metadata_split_groups)) {
    subset_top10 = get_top10_taxa(rare_split_prokaryot_genus[[i]])
    plot_name = paste0('Relative abundance of top 10 prokaryotic genera in ',metadata_split_groups[i],' samples')
    p=make_abundance_barplot(rare_split_prokaryot_genus[[i]],subset_top10,plot_name)
    print(p)
    filename = paste0(output_dir,"/Fig_7-",i,"_abundance_prokaryot_individualtop10_",metadata_split_groups[i],".png")
    export(p,file=filename)
  }
  
  ### Create Abundance barplot for each subset of samples with individual top 10 list, eukaryot ###
  for (i in 1:length(metadata_split_groups)) {
    subset_top10 = get_top10_taxa(rare_split_eukaryot_genus[[i]])
    plot_name = paste0('Relative abundance of top 10 eukaryotic genera in ',metadata_split_groups[i],' samples')
    p=make_abundance_barplot(rare_split_eukaryot_genus[[i]],subset_top10,plot_name)
    print(p)
    filename = paste0(output_dir,"/Fig_8-",i,"_abundance_eukaryot_individualtop10_",metadata_split_groups[i],".png")
    export(p,file=filename)
  }
  
  ### Create Abundance barplot for each subset of samples with individual top 10 list, fungi ###
  for (i in 1:length(metadata_split_groups)) {
    subset_top10 = get_top10_taxa(rare_split_fungi_genus[[i]])
    plot_name = paste0('Relative abundance of top 10 fungal genera in ',metadata_split_groups[i],' samples')
    p=make_abundance_barplot(rare_split_fungi_genus[[i]],subset_top10,plot_name)
    print(p)
    filename = paste0(output_dir,"/Fig_9-",i,"_abundance_fungi_individualtop10_",metadata_split_groups[i],".png")
    export(p,file=filename)
  }
  top_taxa_heatmap(rare_prokaryot,get_top_n_taxa(rare_prokaryot,25),paste0(output_dir,'/prokaryot_heatmap'))
  top_taxa_heatmap(rare_eukaryot,get_top_n_taxa(rare_eukaryot,25),paste0(output_dir,'/eukaryot_heatmap'))
  top_taxa_heatmap(rare_fungi,get_top_n_taxa(rare_fungi,25),paste0(output_dir,'/fungi_heatmap'))
  
  ### Write excel file for kruskal-wallis tests performed on abundance of each genus ###
  genus_comparison_prokaryot = genus_comparison(rare_prokaryot_genus,rare_split_prokaryot_genus)
  genus_comparison_eukaryot = genus_comparison(rare_eukaryot_genus,rare_split_eukaryot_genus)
  genus_comparison_fungi = genus_comparison(rare_fungi_genus,rare_split_fungi_genus)
  filename = paste0(output_dir,"/Table_2_prokaryot_genus_abundances.xlsx")
  
  wb <- createWorkbook("THEJ")
  addWorksheet(wb, "Prokaryot", gridLines = FALSE)
  addWorksheet(wb, "Eukaryot", gridLines = FALSE)
  addWorksheet(wb, "Fungi", gridLines = FALSE)
  
  writeData(wb, sheet = 1, genus_comparison_prokaryot,rowNames = FALSE)
  writeData(wb, sheet = 2, genus_comparison_eukaryot,rowNames = FALSE)
  writeData(wb, sheet = 3, genus_comparison_fungi,rowNames = FALSE)
  saveWorkbook(wb, filename, overwrite = TRUE)
  
  
  
  ### Summary file ###
  summary_file = paste0(output_dir,"/summary_file.txt")
  
  ### Rarefaction stats ###
  threshold_prokaryot = summary_list[[1]]
  threshold_eukaryot = summary_list[[2]]
  threshold_fungi = summary_list[[3]]
  start_nsamples = summary_list[[4]]
  start_ntaxa = summary_list[[5]]
  write(paste0("BION summary\n\n",
               start_nsamples," samples split between ",length(metadata_split_groups)," groups\n"
  ),file = summary_file)
  
  write(paste0("Rarefaction stats. Number of samples included in analysis of each domain"),file = summary_file, append = TRUE)
  rare_matrix = rbind(c("Before rarefaction","-",start_nsamples,start_ntaxa),
                      c("Prokaryot",threshold_prokaryot,nsamples(rare_prokaryot),ntaxa(rare_prokaryot)),
                      c("Eukaryot",threshold_eukaryot,nsamples(rare_eukaryot),ntaxa(rare_eukaryot)),
                      c("Fungi",threshold_fungi,nsamples(rare_fungi),ntaxa(rare_fungi)))
  colnames(rare_matrix) = c("Domain","Rarefaction threshold","Samples included","OTUs included")
  suppressWarnings(write.table(rare_matrix,file=summary_file,append = TRUE, quote = FALSE, sep = "\t",row.names = FALSE))
  
  write(paste0("\n\n\nAlphadiveristy comparison.\nP-values based on Kruskal-Wallis tests.\n"),file = summary_file, append = TRUE)  
  alpha_matrix = cbind(c("Prokaryot","Eukaryot","Fungi"),
                       c(alphadiversity_prokaryot[[3]],alphadiversity_eukaryot[[3]],alphadiversity_fungi[[3]]),
                       c(alphadiversity_prokaryot[[4]],alphadiversity_eukaryot[[4]],alphadiversity_fungi[[4]]))
  colnames(alpha_matrix) = c("Domain","Observed species","Shannon diversity index")
  suppressWarnings(write.table(alpha_matrix,file=summary_file,append = TRUE, quote = FALSE, sep = "\t",row.names = FALSE))
  
  write(paste0("\n\n\nBetadiveristy comparison.\nAnalysis of similarities (ANOSIM) on all groups.\n"),file = summary_file, append = TRUE)
  ANOSIM_matrix = cbind(c("","Prokaryot","Eukaryot","Fungi"),
                        c("R statistic",PCoA_prokaryot[[2]],PCoA_eukaryot[[2]],PCoA_fungi[[2]]),
                        c("p-value",PCoA_prokaryot[[3]],PCoA_eukaryot[[3]],PCoA_fungi[[3]]))
  suppressWarnings(write.table(ANOSIM_matrix,file=summary_file,append = TRUE, quote = FALSE, sep = "\t",row.names = FALSE,col.names =  FALSE))
  
  summary_excel_file = paste0(output_dir,"/summary_tables.xlsx")
  wb <- createWorkbook("THEJ")
  addWorksheet(wb, "Rarefaction stats", gridLines = FALSE)
  addWorksheet(wb, "Alphadiversity p-values", gridLines = FALSE)
  addWorksheet(wb, "ANOSIM stats", gridLines = FALSE)
  
  writeData(wb, sheet = 1, rare_matrix,rowNames = FALSE)
  writeData(wb, sheet = 2, alpha_matrix,rowNames = FALSE)
  writeData(wb, sheet = 3, ANOSIM_matrix,rowNames = FALSE)
  saveWorkbook(wb, summary_excel_file, overwrite = TRUE)
}
