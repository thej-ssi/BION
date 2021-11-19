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

# if(! require("DESeq2")) {
#   source('http://bioconductor.org/biocLite.R')
#   biocLite("DESeq2")
# }
# library('DESeq2')

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

compare_datasets <- function(tsv_input,metadata) {
  BION_IDs = colnames(tsv_input)[8:ncol(tsv_input)]
  metadata_IDs = rownames(metadata)
  missing_BION_IDs = BION_IDs[which(!BION_IDs %in% metadata_IDs)]
  missing_metadata_IDs = metadata_IDs[which(!metadata_IDs %in% BION_IDs)]
  match_IDs = BION_IDs[which(BION_IDs %in% metadata_IDs)]
  print(paste0(length(match_IDs)," sample names are found in both BION data and metadata"))
  print(paste0(length(missing_metadata_IDs)," sample names from BION output are missing from metadata:"))
  if (length(missing_metadata_IDs)>0) {print(missing_metadata_IDs)}
  print(paste0(length(missing_BION_IDs)," sample names from metadata are missing from BION output:"))
  if (length(missing_BION_IDs)>0) {print(missing_BION_IDs)}
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

prune_by_variable <- function(po,variable_name,variable_value) {
  return_po = prune_samples(sample_names(po)[which(get_variable(po,variable_name) %in% variable_value)],po)
}

check_counts <- function(po) {
  print("Number of sequences found in samples")
  print(sort(colSums(otu_table(po))))
  plot(sort(colSums(otu_table(po))),main = "Number of sequences found in samples")
}

make_alphadiversity_object <- function(po,variable_name,plot_title,color_list) {
  groups = levels(factor(get_variable(po,variable_name)))
  if (length(color_list) == length(groups)) {
    col_vec = color_list
  } else {
    col_vec = RColorBrewer::brewer.pal(length(groups),"Set1")
    print(paste0('Number of colors given (', length(color_list) , ') does not match number of levels in variable (', length(groups),')'))
  }
  r <- data.frame(ID=sample_names(po), type=factor(get_variable(po,variable_name)), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon","InvSimpson")))
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
  p3 <- plot_ly(r, y = ~InvSimpson, color = ~type, type = "box", boxpoints = "all", pointpos = -1.5,colors = col_vec) %>%
    layout(title = plot_title,
           #xaxis=list(tickangle = -45),
           yaxis=list(title='Inverse Simpson index'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  shannon_matrix = matrix(ncol = length(groups), nrow = length(groups))
  observed_matrix = matrix(ncol = length(groups), nrow = length(groups))
  simpson_matrix = matrix(ncol = length(groups), nrow = length(groups))
  for (n1 in 1:(length(groups)-1)) {
    print(n1)
    for (n2 in (n1+1):length(groups)) {
      group1 = groups[n1]
      group2 = groups[n2]
      ### Shannon
      vec1 = r$Shannon[which(r$type==group1)]
      vec2 = r$Shannon[which(r$type==group2)]
      wilcox_shannon = wilcox.test(vec1,vec2)
      shannon_matrix[n1,n2] <- wilcox_shannon$p.value
      shannon_matrix[n2,n1] <- wilcox_shannon$p.value
      ### Observed
      vec1 = r$Observed[which(r$type==group1)]
      vec2 = r$Observed[which(r$type==group2)]
      wilcox_observed = wilcox.test(vec1,vec2)
      observed_matrix[n1,n2] <- wilcox_observed$p.value
      observed_matrix[n2,n1] <- wilcox_observed$p.value
      ### Simpson
      vec1 = r$InvSimpson[which(r$type==group1)]
      vec2 = r$InvSimpson[which(r$type==group2)]
      wilcox_simpson = wilcox.test(vec1,vec2)
      simpson_matrix[n1,n2] <- wilcox_simpson$p.value
      simpson_matrix[n2,n1] <- wilcox_simpson$p.value
    }
    
  }
  kruskal_Observed = kruskal.test(r$Observed,r$type)
  kruskal_Shannon = kruskal.test(r$Shannon,r$type)
  kruskal_Simpson = kruskal.test(r$InvSimpson,r$type)
  p_df_shannon = as.data.frame(shannon_matrix)
  p_df_observed = as.data.frame(observed_matrix)
  p_df_simpson = as.data.frame(simpson_matrix)
  rownames(p_df_shannon) = groups
  colnames(p_df_shannon) = groups
  rownames(p_df_observed) = groups
  colnames(p_df_observed) = groups
  rownames(p_df_simpson) = groups
  colnames(p_df_simpson) = groups
  return_list = list("Observed_plot"=p1,"Shannon_plot"=p2,"Simpson_plot"=p3,
                     "Observed_kruskal"=kruskal_Observed$p.value,"Shannon_kruskal"=kruskal_Shannon$p.value,"Simpson_kruskal"=kruskal_Simpson$p.value,
                     "Observed_MWU_mat"=p_df_observed,"Shannon_MWU_mat"=p_df_shannon,"Simpson_MWU_mat"=p_df_simpson)
  print(col_vec)
  return(return_list)
}

make_alphadiversity_object_ggplot <- function(po,variable_name,plot_title,color_list) {
  groups = levels(factor(get_variable(po,variable_name)))
  if (length(color_list) == length(groups)) {
    col_vec = color_list
  } else {
    col_vec = RColorBrewer::brewer.pal(length(groups),"Set1")
    print(paste0('Number of colors given (', length(color_list) , ') does not match number of levels in variable (', length(groups),')'))
  }
  r <- data.frame(ID=sample_names(po), type=factor(get_variable(po,variable_name)), richness=colSums(otu_table(po) > 0), estimate_richness(po,measures = c("Observed","Shannon","InvSimpson")))
  p1 <- ggplot(r,aes(x = type, y = Observed, fill = type)) + geom_boxplot() + geom_point(aes(x=type, y=Observed), position = position_jitter(w = 0.15, h = 0)) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Number of species") + ggtitle(plot_title)
  p2 <- ggplot(r,aes(x = type, y = Shannon, fill = type)) + geom_boxplot() + geom_point(aes(x=type, y=Shannon), position = position_jitter(w = 0.15, h = 0)) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Shannon diversity index") + ggtitle(plot_title)
  p3 <- ggplot(r,aes(x = type, y = InvSimpson, fill = type)) + geom_boxplot() + geom_point(aes(x=type, y=InvSimpson), position = position_jitter(w = 0.15, h = 0)) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Inverse simpson diversity index") + ggtitle(plot_title)
  shannon_matrix = matrix(ncol = length(groups), nrow = length(groups))
  observed_matrix = matrix(ncol = length(groups), nrow = length(groups))
  simpson_matrix = matrix(ncol = length(groups), nrow = length(groups))
  for (n1 in 1:(length(groups)-1)) {
    print(n1)
    for (n2 in (n1+1):length(groups)) {
      group1 = groups[n1]
      group2 = groups[n2]
      ### Shannon
      vec1 = r$Shannon[which(r$type==group1)]
      vec2 = r$Shannon[which(r$type==group2)]
      wilcox_shannon = wilcox.test(vec1,vec2)
      shannon_matrix[n1,n2] <- wilcox_shannon$p.value
      shannon_matrix[n2,n1] <- wilcox_shannon$p.value
      ### Observed
      vec1 = r$Observed[which(r$type==group1)]
      vec2 = r$Observed[which(r$type==group2)]
      wilcox_observed = wilcox.test(vec1,vec2)
      observed_matrix[n1,n2] <- wilcox_observed$p.value
      observed_matrix[n2,n1] <- wilcox_observed$p.value
      ### Simpson
      vec1 = r$InvSimpson[which(r$type==group1)]
      vec2 = r$InvSimpson[which(r$type==group2)]
      wilcox_simpson = wilcox.test(vec1,vec2)
      simpson_matrix[n1,n2] <- wilcox_simpson$p.value
      simpson_matrix[n2,n1] <- wilcox_simpson$p.value
    }
    
  }
  kruskal_Observed = kruskal.test(r$Observed,r$type)
  kruskal_Shannon = kruskal.test(r$Shannon,r$type)
  kruskal_Simpson = kruskal.test(r$InvSimpson,r$type)
  p_df_shannon = as.data.frame(shannon_matrix)
  p_df_observed = as.data.frame(observed_matrix)
  p_df_simpson = as.data.frame(simpson_matrix)
  rownames(p_df_shannon) = groups
  colnames(p_df_shannon) = groups
  rownames(p_df_observed) = groups
  colnames(p_df_observed) = groups
  rownames(p_df_simpson) = groups
  colnames(p_df_simpson) = groups
  return_list = list("Observed_plot"=p1,"Shannon_plot"=p2,"Simpson_plot"=p3,
                     "Observed_kruskal"=kruskal_Observed$p.value,"Shannon_kruskal"=kruskal_Shannon$p.value,"Simpson_kruskal"=kruskal_Simpson$p.value,
                     "Observed_MWU_mat"=p_df_observed,"Shannon_MWU_mat"=p_df_shannon,"Simpson_MWU_mat"=p_df_simpson)
  print(col_vec)
  return(return_list)
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
      print(dim(t(otu_table(po)[,which(!is.na(phen_vec))])))
      dist_obj = vegdist(t(otu_table(po)[,which(!is.na(phen_vec))]),method = "jaccard",binary = TRUE)
      print(dist_obj)
      anosim_test = anosim(dist_obj,grouping = phen_factor,permutations = 1000, distance = "jaccard")
    } else {
      anosim_test = anosim(t(otu_table(po)[,which(!is.na(phen_vec))]),grouping = phen_factor,permutations = 1000, distance = dist_method)
    }
    
    returnlist = list("Plot"=p,"Anosim_results"=anosim_test,"Anosim_R"=anosim_test$statistic,"Anosim_p"=anosim_test$signif)
  } else {
    returnlist = list(p)
  }
  return(returnlist)	
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
make_abundance_barplot_ggplot <- function(po,taxa,plot_title="Relative abundance") {
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
    p = ggplot(melt_top10,aes(x=ID,y=percent,fill=genus)) + geom_bar(position="stack", stat="identity") + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=90)) + ylab("Relative abundance in percent") + ggtitle(plot_title) + 
      scale_fill_manual(values = rev(c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf')))
    return(p)
  }
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
  return_df = as.data.frame(cbind(rep(1:nrow(p_mat)),p_mat,group_means),stringsAsFactors=F)
  #apply(return_df[,9:ncol(return_df)],2, function(e) as.numeric(e))
  colnames(return_df) = c("Rownumber",colnames(tax),mwu_headers,mean_headers)
  rownames(return_df) = rownames(tax)
  return(return_df)
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

raw_read_comparison <- function(po,variable_to_compare) {
  read_sums = as.numeric(colSums(otu_table(po)))
  variable_factor = get_variable(po,variable_to_compare)
  if (!class(variable_factor) == "factor") {
    variable_factor = factor(variable_factor)
  }
  plot_data_frame = data.frame("group" = variable_factor, "read.count" = read_sums)
  print(plot_data_frame)
  p = ggplot()
  groups = levels(variable_factor)
  p_mat = matrix(ncol = length(groups), nrow = length(groups))
  for (n1 in 1:(length(groups)-1)) {
    for (n2 in (n1+1):length(groups)) {
      group1 = groups[n1]
      group2 = groups[n2]
      vec1 = read_sums[which(variable_factor==group1)]
      vec2 = read_sums[which(variable_factor==group2)]
      wilcox_test = wilcox.test(vec1,vec2)
      p_mat[n1,n2] <- wilcox_test$p.value
      p_mat[n2,n1] <- wilcox_test$p.value
    }
    
  }
  colnames(p_mat) = groups
  rownames(p_mat) = groups
  p_df = as.data.frame(p_mat)
  return(list("df"=plot_data_frame,"p.mat"=p_mat))
}

make_heatmap_object <- function(po,top_10_taxa,variable_name,plot_title="Heatmap",color_list=c(),order_by="variable") {
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
    top_10_taxa_ordered = top_10_taxa[order(newnames)]
    top_10_matrix = rbind(otu_table(po)[top_10_taxa_ordered,],other_sum)
    rownames(top_10_matrix) = c(newnames[order(newnames)],"Other")
  } else {
    top_10_matrix = otu_table(po)[top_10_taxa,]
    rownames(top_10_matrix) = newnames
  }
  if (order_by == "clustering") {
    return_object_1 = heatmap.2(top_10_matrix,
                                distfun = vegdist,
                                hclustfun = function(x) hclust(x, method = "ward.D"),
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
  else if (order_by == "none")  {
    return_object = heatmap.2(top_10_matrix,
                                Rowv = FALSE,
                                Colv = FALSE,
                                dendrogram = "none",
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
  else if (order_by == "variable")  {
    top_10_matrix = top_10_matrix[,order(group_color_vector)]
    return_object = heatmap.2(top_10_matrix,
                              Rowv = FALSE,
                              Colv = FALSE,
                              dendrogram = "none",
                              trace = "none",
                              scale = "row",
                              col = colorspace::diverge_hsv(50),
                              margins = c(5,15),
                              ColSideColors = group_color_vector[order(group_color_vector)],
                              #main = "Heatmap showing relative abundance of top 10 genera across all samples",
                              key.title = "",
                              #lwid = 2,
                              labCol = F)
  }
  return(return_object)
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

setup_outdir <- function(dir_name) {
  n = 1
  flag = 0
  while (flag == 0) {
    outdir = paste0(dir_name,'_',n)
    if (file.exists(outdir)) {
      n = n+1
    } else {
      flag = 1
      dir.create(outdir)
    }
  }
  return(outdir)
}

run_cross_sectional_analysis <- function(po, variable_name, color_list, output_folder) {
  output_dir = setup_outdir(output_folder)
  print(paste0("Printing all outputs to ",output_dir))
  sample_data(po)$Group = factor(as.character(as.vector(get_variable(po,variable_name))))
  color_vector = setup_color_vector(po,"Group",color_list)
  
  print("Calculating alphadiversity and printing plots")
  Alphadiv_plot = make_alphadiversity_object_ggplot(po,variable_name = "Group",plot_title = paste0("Alpha diversity grouped by ",variable_name),color_vector)
  filename = paste0(output_dir,"/Fig_1-1_alphadiversity_observed.pdf")
  ggsave(filename = filename,plot = Alphadiv_plot$Observed_plot,device = "pdf",dpi = 300)
  #export(Alphadiv_plot$Observed_plot,file=filename)
  filename = paste0(output_dir,"/Fig_1-2_alphadiversity_shannon.pdf")
  ggsave(filename = filename,plot = Alphadiv_plot$Shannon_plot,device = "pdf",dpi = 300)
  #export(Alphadiv_plot$Shannon_plot,file=filename)
  filename = paste0(output_dir,"/Fig_1-3_alphadiversity_simpson.pdf")
  ggsave(filename = filename,plot = Alphadiv_plot$Simpson_plot,device = "pdf",dpi = 300)
  #export(Alphadiv_plot$Simpson_plot,file=filename)
  
  print("Calculating PCoA and printing plots")
  PCoA_BC = make_PCoA_object(po,variable_name = "Group",plot_title = paste0("PCoA based on Bray Curtis dissimilarity grouped by ",variable_name),color_vector)
  ggsave(filename = paste0(output_dir,"/Fig_2-1_PCoA_BrayCurtis.pdf"),plot = PCoA_BC[[1]],device = "pdf",dpi = 300)
  #PCoA_binary = make_PCoA_object(po,variable_name = "Group",paste0("PCoA based on Binary Jaccard distance grouped by ",variable_name),color_vector,dist_method = "binary")
  #ggsave(filename = paste0(output_dir,"/Fig_2-2_PCoA_binary.pdf"),plot = PCoA_binary[[1]],device = "pdf")
  
  
  ### Barplots ###
  print("Printing barplots")
  po_genus = tax_glom(po,"Genus")
  top10_taxa = get_top_n_taxa(po_genus,10)
  bar_plot = make_abundance_barplot_ggplot(po_genus,taxa = top10_taxa, "Top 10 most abundant genera")
  filename = paste0(output_dir,"/Fig_3-1_barplot_all.pdf")
  ggsave(filename = filename,plot = bar_plot,device = "pdf",dpi = 300)
  #export(bar_plot,file=filename)
  for (i in 1:length(levels(get_variable(po_genus,"Group")))) {
    n = i+1
    var = levels(get_variable(po_genus,"Group"))[i]
    po_sub = prune_by_variable(po_genus,"Group",var)
    var_name = gsub("/","_",var)
    bar_plot = make_abundance_barplot_ggplot(po_sub,taxa = top10_taxa, paste0("Relative abundance of top species in ",var," samples"))
    filename = paste0(output_dir,"/Fig_3-",n,"_barplot_",var_name,".pdf")
    #export(bar_plot,file=filename)
    ggsave(filename = filename,plot = bar_plot,device = "pdf",dpi = 300)
  }
  
  #make_barplot_plus_object(po,top10_taxa,"Group","Top 10 most abundant species",color_vector)
  
  print("Printing heatmap of top 30 genera")
  #Heatmap = make_heatmap_object(po_genus,get_top_n_taxa(po_genus,30),"Group",paste0("Heatmap showing over and underrepresentation of top 30 species, ",variable_name),color_vector)
  pdf(file = paste0(output_dir,"/Fig_4-1_Heatmap.pdf"))
  make_heatmap_object(po_genus,get_top_n_taxa(po_genus,30),"Group",paste0("Heatmap showing over and underrepresentation of top 30 species, ",variable_name),color_vector)
  dev.off()
  
  
  taxa_comparison_df = make_taxa_comparison_object(po_genus,"Group","bonferroni")
  
  pvalue_list = list("Observed_richness"=Alphadiv_plot$Observed_MWU_mat,"Shannon_diversity"=Alphadiv_plot$Shannon_MWU_mat,"Simpson_diversity"=Alphadiv_plot$Simpson_MWU_mat,
                     "PCoA"=data.frame("Type"=c("Bray curtis"),"R.value"=c(PCoA_BC$Anosim_R),"p.value"=c(PCoA_BC$Anosim_p)),
                     "Genus_abundance_comparison"=taxa_comparison_df)
  print("Printing p values to excel sheets")
  write.xlsx(pvalue_list,paste0(output_dir,"/p_value_tables.xlsx"),row.names=TRUE)
  
  #return(list(Alphadiv_plot,PCoA_BC,PCoA_binary,bar_plot,Heatmap,taxa_comparison_df))
}


