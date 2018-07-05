###########################
#####                 #####
#####################  Set variables  #####################
#####                 #####
###########################

##### BION table and metadata files #####
setwd("/Volumes/data/MPV/THEJ/BION/Sarcoidosis/") ### NEEDED INPUT. Folder containing input files
input_file = "counts_2.txt"                    ### NEEDED INPUT. Name of file containing sequence counts (.txt og .xlsx)
metadata_file = "metadata_fixed.xlsx"                           ### NEEDED INPUT. Name of metadata file (.txt og .xlsx)

##### metadata variable for sub set plots #####
metadata_split_variable = "Diagnose"                         ### NEEDED INPUT. Header of the column that indicates which group the sample is in ###
metadata_split_groups = c("Sarcoidosis","Non-specific inflammation","Cancer")           ### NEEDED INPUT. Names of the different groups to be compared



### Load source file ###
#source("P:/MPV/THEJ/BION/FABIAN_backup/automate_functions.R")
#source("https://raw.githubusercontent.com/thej-ssi/BION/master/automate_functions.R")

##### Load input files #####
tsv_input = load_data(input_file)
metadata = load_metadata(metadata_file,metadata_split_variable)

### Create phyloseq objects for prokaryot, species, eukaryot species, and funal species ###
phylo_objects<-setup_phylo_object(tsv_input,metadata)
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
prokaryot_rarefaction_threshold = 970         ### NEEDED INPUT, rarefaction threshold for prokaryotes
eukaryot_rarefaction_threshold = 452           ### NEEDED INPUT, rarefaction threshold for eukaryotes
fungi_rarefaction_threshold = 445              ### NEEDED INPUT, rarefaction threshold for fungi


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


po = rare_fungi
top_10_taxa = get_top_n_taxa(po,30)

groups = levels(sample_data(po)$sample_groups)
group_colors = RColorBrewer::brewer.pal(n=3,name="Set1")
color_names = c("Red","Blue","Green","Purple","Orange","Yellow","Brown","Pink","Grey")
group_color_vector = as.vector(sample_data(po)$sample_groups)

color_table = matrix(ncol=2,nrow=0)
for (i in 1:length(groups)) {
  group_color_vector[group_color_vector==groups[i]] = group_colors[i]
  color_table= rbind(color_table,c(groups[i],color_names[i]))
}

image_filename = paste0(output_dir,"/Fig_10a_phylum_heatmap_top10_genus.png")
legend_filename = paste0(output_dir,"/Fig_10_heatmap_color_legend.txt")
write.table(color_table,file = legend_filename,sep="\t",quote = F,col.names = F,row.names = F)

other_sum = colSums(otu_table(po)[-top_10_taxa,])
top_10_matrix = rbind(otu_table(po)[top_10_taxa,],other_sum)
top_10_names = c(paste0(tax_table(po)[top_10_taxa,6]," ",tax_table(po)[top_10_taxa,7]),"Other")
rownames(top_10_matrix)=top_10_names

#png(filename = image_filename, width = 1000, height = 700, res=72)
heatmap.2(top_10_matrix,
          trace = "none",
          scale = "row",
          col = colorspace::diverge_hsv(50),
          margins = c(5,22),
          ColSideColors = group_color_vector,
          #main = "Heatmap showing relative abundance of top 10 genera across all samples",
          key.title = "",
          #lwid = 2,
          cexRow = 1.8,
          labCol = F)
#dev.off()

top_10_test = get_top10_taxa(rare_prokaryot)

top_taxa_heatmap(rare_prokaryot,top_10_test)


rare_prokaryot_phylum = tax_glom(rare_prokaryot,taxrank = "Phylum")

top10_phyla = get_top10_taxa(rare_prokaryot_phylum)

top_taxa_heatmap(rare_prokaryot_phylum,top10_phyla,"heatmap_prokaryot_phylum")

top10_species = get_top10_taxa(rare_prokaryot)
make_abundance_barplot_2(rare_prokaryot,top10_species,plot_name = 'test')

make_abundance_barplot_2(rare_prokaryot_genus,get_top10_taxa(tax_glom(rare_prokaryot,taxrank = "Genus")),plot_name = 'test')

top20_species = get_top_n_taxa(rare_prokaryot,20)
top_taxa_heatmap(rare_prokaryot,top20_species,"heatmap_prokaryot_species")


top_taxa_heatmap(rare_fungi_genus,get_top_n_taxa(rare_fungi_genus,10),"heatmap_fungi_species_genus")

top_taxa_heatmap(rare_prokaryot,get_top_n_taxa(rare_prokaryot,25),"heatmap_prokaryot_species")


top_taxa_heatmap(rare_prokaryot_genus,get_top_n_taxa(rare_prokaryot_genus,25),"heatmap_prokaryot_genus")


top_taxa_heatmap(rare_prokaryot_genus,get_top_n_taxa(rare_prokaryot_genus,10),"heatmap_prokaryot_genus")

po = rare_fungi_genus
top_10_taxa = get_top_n_taxa(po,10)
plot_name = "Relative abundance of top 10 fungal genera"

make_abundance_barplot_with_heatmap_tiles <- function(po,top_10_taxa,plot_name) {
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
    
    heatmap_matrix = as.matrix(dd[1:2,fit$order])
    groups = levels(sample_data(po)$sample_groups)
    group_colors = RColorBrewer::brewer.pal(n=length(groups),name="Set1")
    color_names = c("Red","Blue","Green","Purple","Orange","Yellow","Brown","Pink","Grey")
    group_color_vector = as.vector(sample_data(po)$sample_groups)
    for (i in 1:length(groups)) {
      group_color_vector[which(group_color_vector==groups[i])] = group_colors[i]
    }
    heatmap.2(heatmap_matrix,
              Colv = F,
              Rowv = F,
              dendrogram = "none",
              ColSideColors = group_color_vector,
              key = F)
    return(p)
  }
}


make_abundance_barplot_with_heatmap_tiles(rare_prokaryot_genus,get_top_n_taxa(rare_prokaryot_genus,10),plot_name = 'Relative abundance of top 10 prokaryot genera')


make_abundance_barplot_with_heatmap_tiles(rare_fungi_genus,get_top_n_taxa(rare_fungi_genus,10),plot_name = 'Relative abundance of top 10 fungal genera')
