library(plotly)
if(! require("RCurl")) {install.packages("RCurl")}
source('http://bioconductor.org/biocLite.R')

if(! require("phyloseq")) {
  biocLite('phyloseq')
}
library(phyloseq)
library(vegan)
library(tibble)
library(reshape)
library(plotly)
library(gplots)
library(yaml)
library(RCurl)
library(backports)
#library(DESeq2)
library(xlsx)
library(webshot)

setwd("H:/SFSDI/MIAFD/Mikrobiom/Runs/2018/MykParLab124/R")

input_file = "16S.species.tsv"

loaded_input = load_data(input_file)

tsv_input = loaded_input[[1]]
read_counts = loaded_input[[2]]
mapped_counts = loaded_input[[3]]

box_plot_sample(as.numeric(read_counts),as.numeric(read_counts[30]),'test')

box_plot_sample(as.numeric(mapped_counts),as.numeric(mapped_counts[30]),'test')


count_vector = as.numeric(mapped_counts)


tax_vector = as.vector(tsv_input$Taxonomic.groups)
otu_table = tsv_input[,!colnames(tsv_input) %in% c("Taxonomic.groups","Row.max")]
tax_table = matrix(nrow = 0, ncol = 7)
for (i in 1:length(tax_vector)) {
  tax_split = strsplit(tax_vector[i],"; ",fixed=TRUE)
  tax_vec = c(tax_split[[1]][1],tax_split[[1]][2],tax_split[[1]][3],tax_split[[1]][4],tax_split[[1]][5],tax_split[[1]][6],tax_split[[1]][7])
  tax_table = rbind(tax_table,tax_vec)
}
rownames(otu_table) = tax_vector
rownames(tax_table) = tax_vector
colnames(tax_table) = c("Domain","Phylum","Class","Order","Family","Genus","Species")

po = phyloseq(otu_table(otu_table,taxa_are_rows=TRUE),tax_table(tax_table))

po_genus = tax_glom(po,taxrank = "Genus")
po_phylum = tax_glom(po,taxrank = "Phylum")

load_data <- function(input_file) {
  
  tsv_input = read.table(input_file,sep = '\t', comment.char = "", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (substr(tsv_input[1,1],1,1) == '#') {
    colcount = ncol(tsv_input)
    sample_vector = colnames(tsv_input)[c(1:(colcount-6))]
    sample_vector[1] = substr(sample_vector[1],4,nchar(sample_vector[1]))
    read_count_vector = tsv_input[2,c(1:(colcount-6))]
    read_count_vector[1] = substr(read_count_vector[1],4,nchar(read_count_vector[1]))
    mapped_count_vector = tsv_input[3,c(1:(colcount-6))]
    mapped_count_vector[1] = substr(mapped_count_vector[1],4,nchar(mapped_count_vector[1]))
    counts = tsv_input[c(4:nrow(tsv_input)),c(1:(colcount-6))]
    tax_vector = as.vector(tsv_input[c(4:nrow(tsv_input)),ncol(tsv_input)])
    colnames(counts) = sample_vector
    names(read_count_vector) = sample_vector
    names(mapped_count_vector) = sample_vector
    
    new_counts = matrix(nrow = nrow(counts), ncol = 0)
    for (col in 1:ncol(counts)) {
      new_vector = as.numeric(as.vector(counts[,col]))
      new_counts = cbind(new_counts,new_vector)
    }
    colnames(new_counts) = sample_vector
    tsv_input = as.data.frame(new_counts)
    tsv_input$Taxonomic.groups = tax_vector
  }
  returnlist = list(tsv_input,read_count_vector,mapped_count_vector)
  return(returnlist)
}

box_plot_sample <- function(count_vector,sample_count,plot_name) {
  p1 <- plot_ly(x = count_vector, type = "box", boxpoints = "all", pointpos = -1.5, name = "Full run") %>%
    add_trace(x = sample_count, name = 'Isolate', type = "scatter") %>%
    layout(title = plot_name,
           #xaxis=list(tickangle = 90),
           #yaxis=list(title='Number of observed OTUs'),
           margin = list(l=50,r=50,b=100,t=50),
           showlegend = FALSE)
  print(p1)
}

setup_phylo_object <- function(tsv_input) {
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
  Aspergillaceae_remove = collapse_Aspergillaceae(tax_table,otu_table)
  tax_table = Aspergillaceae_remove[[1]]
  otu_table = Aspergillaceae_remove[[2]]
  tax_prokaryot = as.matrix(tax_table[(tax_table[,1]=="d__Bacteria"|tax_table[,1]=="d__Archaea") & !tax_table[,2]=="p__Cyanobacteria/Chloroplast",])
  otu_prokaryot = data.matrix(otu_table[(tax_table[,1]=="d__Bacteria"|tax_table[,1]=="d__Archaea") & !tax_table[,2]=="p__Cyanobacteria/Chloroplast",])
  tax_eukaryot = as.matrix(tax_table[tax_table[,1]=="d__Eukaryota" & !tax_table[,3]=="c__Mammalia",])
  otu_eukaryot = data.matrix(otu_table[tax_table[,1]=="d__Eukaryota" & !tax_table[,3]=="c__Mammalia",])
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

  if (nrow(tax_prokaryot) > 0) {
    rownames(tax_prokaryot) <- paste0("OTU", 1:nrow(tax_prokaryot))
    rownames(otu_prokaryot) = rownames(tax_prokaryot)
    colnames(tax_prokaryot) = tax_ranks
    po_prokaryot = phyloseq(tax_table(tax_prokaryot),otu_table(otu_prokaryot,taxa_are_rows = TRUE))
  } else {
    po_prokaryot = 0
  }
  if (nrow(tax_eukaryot) > 0) {
    rownames(tax_eukaryot) <- paste0("OTU", 1:nrow(tax_eukaryot))
    rownames(otu_eukaryot) = rownames(tax_eukaryot)
    colnames(tax_eukaryot) = tax_ranks
    po_eukaryot = phyloseq(tax_table(tax_eukaryot),otu_table(otu_eukaryot,taxa_are_rows = TRUE))
  } else {
    po_eukaryot = 0
  }
  if (length(fungi_rows)>0) {
    tax_fungi = tax_table[fungi_rows,]
    otu_fungi = otu_table[fungi_rows,]
    rownames(tax_fungi) <- paste0("OTU", 1:nrow(tax_fungi))
    rownames(otu_fungi) = rownames(tax_fungi)
    colnames(tax_fungi) = tax_ranks
    po_fungi = phyloseq(tax_table(tax_fungi),otu_table(otu_fungi,taxa_are_rows = TRUE))
  } else {
    po_fungi = 0
  }
  returnlist = list(po_prokaryot,po_eukaryot,po_fungi)
  return(returnlist)
}

collapse_Aspergillaceae <- function(tax_table,otu_table) {
  Aspergillaceae_rows = which(tax_table[,5]=="f__Aspergillaceae")
  Trichocomaceae_rows = which(tax_table[,5]=="f__Trichocomaceae")
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
        tax_table[i,]=c(as.vector(tax_table[i,1:4]),"f__Trichocomaceae",as.vector(tax_table[i,6:7]))
      }
      
    }
    if (length(remove_rows)>0){
      tax_table = tax_table[-remove_rows,]
      otu_table = otu_table[-remove_rows,]
    }
  }
  Aspergillaceae_rows = which(tax_table[,5]=="f__Aspergillaceae")
  for (i in Aspergillaceae_rows) {
    asper_tax = as.vector(tax_table[i,])
    asper_tax[5] = 'f__Trichocomaceae'
    tax_table[i,] = asper_tax
    print(asper_tax)
  }
  return_list = list(tax_table,otu_table)
}

make_genus_piechart <- function(po,sample,top10_taxa,plot_name) {
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
        newname = paste0(substr(tax_vector[6],4,nchar(tax_vector[6])),' ',substr(tax_vector[7],4,nchar(tax_vector[7])))
      } else if (!tax_vector[6]=="g__unclassified") {
        newname = substr(tax_vector[6],4,nchar(tax_vector[6]))
      } else if (!tax_vector[5]=="f__unclassified") {
        newname = substr(tax_vector[5],4,nchar(tax_vector[5]))
      } else if (!tax_vector[4]=="o__unclassified") {
        newname = substr(tax_vector[4],4,nchar(tax_vector[4]))
      } else if (!tax_vector[3]=="c__unclassified") {
        newname = substr(tax_vector[3],4,nchar(tax_vector[3]))
      } else {
        newname = substr(tax_vector[2],4,nchar(tax_vector[2]))
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

otu_table = otu_table(po_phylum)
tax_table = tax_table(po_phylum)

isolate_counts = otu[,3]
otus_sorted<-sort.list(isolate_counts, decreasing=T)
top9 = otus_sorted[1:9]
top9 = top9[which(isolate_counts[top9]>0)]
other = otus_sorted[-top9]
isolate_df = cbind(c(as.vector(tax_table[top9,2]),'Other'),c(as.vector(isolate_counts[top9]),sum(as.vector(isolate_counts[other]))))


p <- plot_ly(as.data.frame(otu_table),values = otu_table[1:10,],labels = tax_table[1:10,6], type = 'pie') %>%
  layout(title = 'United States Personal Expenditures by Categories in 1960',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
p


USPersonalExpenditure <- data.frame("Categorie"=rownames(USPersonalExpenditure), USPersonalExpenditure)
data <- USPersonalExpenditure[,c('Categorie', 'X1960')]

test_df = data.frame(phylum=c('Actinobacteria','Bacteroidetes','Euryarchaeota','Firmicutes','Synergistetes'), counts=c(5,50,10,45,10))

p <- plot_ly(test_df, labels = ~phylum, values = ~counts, type = 'pie',
             textposition = 'outside',
             marker = list(showlegend = TRUE)
             insidetextfont = list(color = '#FFFFFF'),
             textinfo = 'label+percent') %>%
  layout(title = 'Bacterial phylum composition of sample 217',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

p
