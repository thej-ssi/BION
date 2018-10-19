###################################################
## Load source script with functions from Github ##
###################################################

source("https://raw.githubusercontent.com/thej-ssi/BION/master/automate_functions.R")

########################
## Read into phyloseq ##
########################

##### Set working directory #####
setwd("/Users/thej/Documents/DFPT/")

##### Load metadata #####
m<-read.table("metadata.txt", sep="\t", h=T, stringsAsFactors=F, row.names=1)

##### Load BION data #####
d<-load_data("16S.species.tsv")

##### Setup phyloseq object #####
po = setup_phylo_object(d,m)


##### This returns a list of 3 objects #####
## po[[1]] contains prokaryots
## po[[2]] contains eukaryots
## po[[3]] contains fungi

### This dataset only has 16s prokaryotic data, so we just use the first one
po_raw = po[[1]]

### This combines some taxa that are identical but listed as seperate in BION output
po_raw = combine_duplicates_from_phylo_object(po_raw)


### Have a look
check_counts(po_raw)

### Rarefy data set to desired threhold
rarefaction_threshold = 26602
po = rarefy_even_depth(po_raw,sample.size = rarefaction_threshold,rngseed = 100)


##### Phyloseq objects have three components in them

### a taxonomy table that can be accessed with tax_table(po)
### an OTU table that can be accessed with otu_table(po)
### a metadata table that can be accessed with sample_data(po)

### have a look at the first line of the metadata:
sample_data(po)[1,]


#### We're going to include an extra variable in the metadata, "D.fragilis", which is set depending on the qPCR measurement in the table

sample_data(po)$D.fragilis = sample_data(po)$D.fragilis.Ct
sample_data(po)$D.fragilis[which(m$D.fragilis.Ct<10)] = "Negative"
sample_data(po)$D.fragilis[which(m$D.fragilis.Ct>=10)] = "Positive"

### Right now our data is split into OTUs at species level, we can also look at other taxonomic levels
po_genus = tax_glom(po,taxrank = "Genus")
po_phylum = tax_glom(po,taxrank = "Phylum")


### We can also set up phyloseq objects for subsets of data, like all controls, all male, etc.

male_samples = rownames(m)[which(m$Gender=="M")]
female_samples = rownames(m)[which(m$Gender=="F")]
all_time_samples = rownames(m)[which(m$Time.check==1)]
treatment_samples = rownames(m)[which(m$Treatment==1)]
control_samples = rownames(m)[which(m$Treatment==0)]

t0_samples = rownames(m)[which(m$Sample.time==0)]
t1_samples = rownames(m)[which(m$Sample.time==1)]
t2_samples = rownames(m)[which(m$Sample.time==2)]
t3_samples = rownames(m)[which(m$Sample.time==3)]
t4_samples = rownames(m)[which(m$Sample.time==4)]
t5_samples = rownames(m)[which(m$Sample.time==5)]
t6_samples = rownames(m)[which(m$Sample.time==6)]

age1_samples = rownames(m)[which(m$Age.group==1)]
age2_samples = rownames(m)[which(m$Age.group==2)]
age3_samples = rownames(m)[which(m$Age.group==3)]

sample_data(po)$Group = sample_data(po)$Treatment
sample_data(po)$Group[which(sample_data(po)$Treatment==1)] = "Treatment"
sample_data(po)$Group[which(sample_data(po)$Treatment==0)] = "Control"

po_treatment = prune_samples(treatment_samples,po)
po_control = prune_samples(control_samples,po)
po_t0 = prune_samples(t0_samples,po)
po_t1 = prune_samples(t1_samples,po)
po_t2 = prune_samples(t2_samples,po)
po_t3 = prune_samples(t3_samples,po)
po_t4 = prune_samples(t4_samples,po)
po_t5 = prune_samples(t5_samples,po)

po_timecheck = prune_samples(all_time_samples,po)



### Before we continue let's have a quick look at some colors
color_vector = RColorBrewer::brewer.pal(9,"Set1")

test_color_tile(color_vector)


### Let's try a simple alphadiversity plot with treatment and control samples, we'll use green and (1st and 3rd in the color vector)

Alphadiv_plot = make_alphadiversity_object(po,variable_name = "Group",plot_name = "Alphadiversity of control and treatment samples",color_vector[c(1,3)])
Alphadiv_plot


#### "Alphadiv_plot" is a list of 5 different components

# 1. A plot with number of observed OTUs
# 2. A plot with Shannon diversity
# 3. p-value from kruskal-wallis test performed on observed OTUs
# 4. p-value from kruskal-wallis test performed on Shannon diversity
# 5. matrix of p-values from pairwise Mann-Whitney U tests performed on Shannon diversity of the different groups (here we only have 2 groups, so a simple matrix)


#### Let's try a PCoA plot with the same color setup

PCoA_plot = make_PCoA_object(po,variable_name = "Group",plot_title = "PCoA plot of control and treatment samples",color_vector[c(1,3)])
PCoA_plot


### Igen får vi en liste med forskellige ting
# 1. er selve PCoA plottet
# 2. er R værdien fra en anosim test (indikerer hvor stor adskillelse der er mellem grupperne)
# 3. er p-værdien fra en anosim test (indikerer hvorvidt adskillelsen er signifikant)


### Lad os se på nogle barplots

# To make barplots we need a list of which taxa to include. To get the most abundant taxa in the data set use
top10_taxa = get_top_n_taxa(po,10)

# And input that into the barplot function
bar_plot = make_abundance_barplot(po,taxa = top10_taxa, "Top 10 most abundant species")
bar_plot


### If you're in for a little paint job we can add some things to the barplot

make_barplot_plus_object(po,top10_taxa,"Group","Top 10 most abundant species",color_vector[c(1,3)])

### This will create two new plots under "Viewer" to the right, and one plot under "Plots"
# The first plot under Viewer (if you go one back) is the same barplot as before
# Second plot under viewer is the alphadiversity of all the samples, listed in the same order as they appear in the plot
# Under "Plots" there will be a heatmat. The heatmap is irrelevant, but the tiles on top are ordered according to the bars and indicate which group each sample belongs to


### Heatmaps

Heatmap = make_heatmap_object(po,get_top_n_taxa(po,30),"Group","Heatmap showing over- and underrepresentation of top 30 species",color_vector[c(1,3)])
Heatmap


### We can also have a look at specific OTUs

taxa_comparison_df = make_taxa_comparison_object(po,"Group","bonferroni")

# this performs a Mann-Whitney U test on all taxa and makes a list
# Lets have a look at just the lines in the dataframe where the multiple correction adjusted p-value is below 0.05

taxa_comparison_df[which(as.vector(taxa_comparison_df$`p adjusted`)<0.05),]

# The number in the rownumber column indicates which row in our otu_table we find that taxa. We can make a boxplot of the relative abundance for that specific OTU
# For instance to have a look at Alistipes putredinis (Rownumber 160) use

OTU_boxplot = make_OTU_boxplot_object(po,OTU = 160,variable_name = "Group", plot_name = "Relative abundance distribution",color_vector[c(1,3)])
OTU_boxplot


