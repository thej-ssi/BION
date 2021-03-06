
##### Check if dependencies are installed and load packages #####
if(! require("RCurl")) {install.packages("RCurl")}
source('http://bioconductor.org/biocLite.R')

if(! require("phyloseq")) {
  biocLite('phyloseq')
}
if(! require("DESeq2")) {
  biocLite("DESeq2")
}
library('phyloseq')
require('vegan')
require('tibble')
require('reshape')
require('plotly')
require('gplots')
require('yaml')
require('RCurl')
require('backports')
require('DESeq2')
require('openxlsx')
require('webshot')
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


### Set up two phyloseq objects, one for prokaryot and one for eukaryot species ###
setup_phylo_object <- function(tsv_input,metadata) {
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
  returnlist = list(po_prokaryot,po_eukaryot,po_fungi)
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

make_abundance_barplot <- function(po,top_10_taxa,plot_name) {
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
  r <- data.frame(ID=sample_names(po), type=factor(get_variable(po)$sample_groups), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon")))
  p1 <- plot_ly(r, y = ~Observed, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5) %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = 90),
           yaxis=list(title='Number of observed OTUs'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  p2 <- plot_ly(r, y = ~Shannon, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5) %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = -45),
           yaxis=list(title='Shannon diversity index'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  kruskal_Observed = kruskal.test(r$Observed,r$type)
  kruskal_Shannon = kruskal.test(r$Shannon,r$type)
  print(p1)
  print(p2)
  return_list = list(p1,p2,kruskal_Observed$p.value,kruskal_Shannon$p.value)
  return(return_list)
}


set_api_key <- function() {
  Sys.setenv("plotly_username"="thej-ssi")
  Sys.setenv("plotly_api_key"="gFKgrgfaQjKs1GanZdA7")
}

make_PCOA_plot <- function(po,plotname) {
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

check_counts_both <- function(po_prokaryot,po_eukaryot) {
  print("Number of prokaryot sequences found in samples")
  print(sort(colSums(otu_table(po_prokaryot))))
  plot(sort(colSums(otu_table(po_prokaryot))),main = "Number of prokaryot sequences found in samples")
  print("Number of eukaryot sequences found in samples")
  print(sort(colSums(otu_table(po_eukaryot))))
  plot(sort(colSums(otu_table(po_eukaryot))),main = "Number of eukaryot sequences found in samples")
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
  
  image_filename = paste0(output_dir,"/",filename,".png")
  legend_filename = paste0(output_dir,"/",filename,"__color_key.txt")
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
