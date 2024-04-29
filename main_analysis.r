# based on https://software.cqls.oregonstate.edu/tips/posts/phyloseq-bug-meeting-presentation-fall-2019/
rm(list=ls(all=TRUE))
library(ggplot2)
library(decontam)
library(phyloseq)
library(Biostrings)
library(DESeq2)
library("ape")
library(knitr)
library(stringr)
library("vegan")
library('ggpubr')
library('cowplot')
library('ampvis2')
library('ggvegan')
library('ggbiplot')
library("dichromat")
library("RColorBrewer")
library('UpSetR')
library(rstatix)
library("pheatmap") 
library(dplyr)
library("vsn")
library(writexl)
library("grid")
library("gridExtra")
library(readxl)
library('biomformat')
mac.path <- "/Users/riccardo/Library/CloudStorage/GoogleDrive-riccardo.piccinno@unipv.it/.shortcut-targets-by-id/1OLHPrObneyH2cOmEvB_ahqndO3QTORXC/"
win.path <- "G:/Il mio Drive/PhD.AES.Dsuz.Riccardo.Piccinno/"
load(paste0(mac.path,"H_halys/metaBUG/Sequenze_MiSeq/cutadapt_dada2.RData"))
setwd(paste0(mac.path,"H_halys/metaBUG/Sequenze_MiSeq"))
source(paste0(mac.path,"H_halys/metaBUG/Sequenze_MiSeq/functions.r"))

###### dataset preparation #####
samples.out <- rownames(seqtab.nochim)
subject <- str_replace_all(samples.out, 'K-', 'K')
subject <- sapply(strsplit(subject, "-"), `[`, 2)
population <- str_replace_all(subject, '[0-9]', '')
population <- str_replace_all(population, 'Kex', 'K')
population <- str_replace_all(population, 'Kpcr', 'K')
n <- str_replace_all(subject, '[A-z]', '')
n[5] <- 1
samdf <- data.frame(Subject=subject, Population=population, Number=n, 
                    quant_reading=as.data.frame(track)$nonchim)
samdf$Overwintered <- "No"
samdf$Overwintered[samdf$Population=='OV'] <- "Yes"
rownames(samdf) <- samples.out
ov.names <- c('TNS','TNS','TNN','SMA','TNS','SMA','TNN',
              'VN','VN','TNN','TNN','VN','VN','SMA',
              'TNS','TNN','VN','TNS','VN','SMA','TNN')
samdf$Population[6:(6+20)] <- ov.names 
samdf$Population[samdf$Population == "TNN"] <- "TN"
samdf$Population[samdf$Population == "VN"] <- "DN"
samdf$Population[samdf$Population == "TNS"] <- "BS"
samdf$Population[samdf$Population == "SMA"] <- "SM"
samdf$Subject <- gsub("^TNS", "BS", samdf$Subject)
samdf$Subject <- gsub("^SMA", "SM", samdf$Subject)
samdf$Subject <- gsub("^TNN", "TN", samdf$Subject)
samdf$Subject <- gsub("^VN", "DN", samdf$Subject)
# Generating a phyloseq object with all phylogenetic information on OTUs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
tree <- rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna, tree)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
df$is.neg <- df$Population == "K"
sample_data(ps)$is.neg <- sample_data(ps)$Population == "K"
ggplot(data=df, aes(x=Index, y=LibrarySize, color=is.neg)) + geom_point()
contamdf <- isContaminant(ps, method="either", conc="quant_reading", neg='is.neg', threshold=0.05)
table(contamdf$contaminant)
which(contamdf$contaminant)
# The default threshold for a contaminant is that it reaches a probability of 
# 0.1 in the statistical test being performed. In the prevalence test there is 
# a special value worth knowing, threshold=0.5, that will identify as contaminants 
# all sequences thare are more prevalent in negative controls than in positive samples.
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Population == "K", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Population != "K", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ps.contam <- prune_taxa(contamdf$contaminant, ps)
ps.noncontam <- prune_taxa(!contamdf$contaminant, ps)
ps.noncontam <- prune_samples(sample_data(ps.noncontam)$Population != "K", ps.noncontam)
ps <- ps.noncontam
ncolors <- length(levels(as.factor(sample_data(ps)$Population)))-1
primary_color_list <- brewer.pal(ncolors, 'Dark2') 
ps <- prune_samples(sample_data(ps)$Population != "VIG", ps)
sample_names(ps) <- sample_data(ps)$Subject
ps.prop <- transform_sample_counts(ps, function(x) x/sum(x)*100)

###### save OTU table ######
OTU1 <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
OTUdf <- data.frame(t(OTU1)) # taxa rows, samples columns
write_xlsx(cbind(" "=rownames(OTUdf), OTUdf),'OTU_table.xlsx')
OTUbiom <- make_biom(OTUdf)
write_biom(OTUbiom, 'OTU_table.biom')

##### Mantel test (Beta diversity test for population connectivity) #####
otu_dist <- vegdist(OTU1, method = "bray") # compute otu distance
loc <- read_excel(paste0(mac.path,"H_halys/metaBUG/Locations.xlsx")) # read locations' coordinates
loc_data <- left_join(sample_data(ps),loc, by = "Population") # add location coordinates for each sample
loc_coord <- as.data.frame(cbind(as.numeric(loc_data$Coordinate1),as.numeric(loc_data$Coordinate2))) # make a dataframe containing only coordinates info
library(geosphere) # load package to compute geographic distances
geo_dist <- distm(loc_coord) # perform geographic distances among points
mantel(otu_dist, geo_dist, method = "pearson", permutations = 9999) # mantel test

##### save OTU sequences ######
ps %>%
  refseq() %>%
  Biostrings::writeXStringSet("./asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

refseq.df <- as.data.frame(refseq(ps))

##### Pantoea numbers #####
ps2.bar <- ps
ps2.bar.merge.comp <- transform_sample_counts(ps2.bar,  
                                              function(x) x/sum(x) * 100)
ps2.bp.Genus <- tax_glom(ps2.bar.merge.comp, 'Genus', NArm = F)
pantoea <- subset_taxa(ps2.bp.Genus, Genus=="Pantoea")
pant.high <- prune_samples(sample_data(pantoea)$Population == "SM", pantoea)
length(otu_table(pant.high))
otu_table(pant.high)>50
pant.high <- prune_samples(sample_data(pantoea)$Population == "BS", pantoea)
length(otu_table(pant.high))
otu_table(pant.high)>50
pant.high <- prune_samples(sample_data(pantoea)$Population == "TN", pantoea)
length(otu_table(pant.high))
otu_table(pant.high)>50
pantoea <- subset_taxa(ps2.bp.Genus, Genus=="Commensalibacter")
pant.high <- prune_samples(sample_data(pantoea)$Population == "DN", pantoea)
length(otu_table(pant.high))
otu_table(pant.high)>50
pant.high <- prune_samples(sample_data(pantoea)$Population == "TN", pantoea)
length(otu_table(pant.high))
otu_table(pant.high)>50
pantoea <- subset_taxa(ps2.bp.Genus, Genus=="Yokenella")
pant.high <- prune_samples(sample_data(pantoea)$Population == "TN", pantoea)
length(otu_table(pant.high))
otu_table(pant.high)>50

# check pathogens
pathogens <- c("Agrobacterium","Burkholderia","Erwinia","Pectobacterium","Pseudomonas","Ralstonia","Xanthomonas")
for (p in pathogens){
  tryCatch({
    path <- subset_taxa(ps2.bp.Genus, Genus==p)
    abd <- sum(otu_table(path))
    n <- length(otu_table(path)[otu_table(path)>0])
    print(paste(p,"abundance:",abd,"samples:",n))
    if (n > 0){
      x <- otu_table(path)[otu_table(path)>0]
      print(rownames(x))
      x <- sample_data(path)[rownames(sample_data(path)) %in% rownames(x),]
      print(paste("TN",dim(x[x=="TN"])[1]))
      print(paste("BS",dim(x[x=="BS"])[1]))
      print(paste("SM",dim(x[x=="SM"])[1]))
      print(paste("DN",dim(x[x=="DN"])[1]))
      }
  }, error = function(e) {
    print(paste(p,0,"abundance:",0,"samples:",0))
  })
  
}

ps_pseudomonas <- subset_taxa(ps, Species == "cichorii")
pseudomonas <- unique(tax_table(ps_pseudomonas)[,c("Genus","Species","Species.1")])

ps.gen <- tax_glom(ps,"Genus",NArm=F)
ps.gen.prop <- transform_sample_counts(ps.gen, function(x) x/sum(x)*100)
filtered_ps <- prune_taxa(taxa_sums(ps.gen.prop) >= 0.3, ps.gen.prop)

my_plot_bar <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                         facet_grid = NULL) {
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack") + xlab('') + labs(x='') + scale_fill_brewer(palette="Pastel1")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 1))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}
clade = "Genus"
ps.prop.renamed <- ps.gen.prop
cabd <- sort(taxa_sums(ps.prop.renamed), decreasing = T)
# Lets remove some low abundance phyla
renamed.clade <- names(cabd[4:length(cabd)])
tax_table(ps.prop.renamed)[renamed.clade, clade] <- "Other"
png(paste('./pop_phylobars_full_',clade,'.png',sep=''), width=9.5,units="in", height=5, res=300)
print(my_plot_bar(ps.prop.renamed, fill=clade, x="Subject") + facet_wrap(~Population, scales="free_x", nrow=1)+guides(x='none')+ theme(text = element_text(size = 16)))
dev.off()

png(paste('./ov_phylobars_full_',clade,'.png',sep=''), width=9.5,units="in", height=5, res=300)
print(my_plot_bar(ps.prop.renamed, fill=clade, x="Subject") + facet_wrap(~Overwintered, scales="free_x", nrow=1)+guides(x='none')+ theme(text = element_text(size = 16)))
dev.off()

####### population analysis ######
pop.div.plots <-alpha_diversity_pop(ps)
pop.div.plots_phy <-alpha_diversity_pop(ps,"Phylum")
pop.div.plots_gen <-alpha_diversity_pop(ps,'Genus')
all.div.plots <-alpha_diversity_ov(ps)
all.div.plots_phy <-alpha_diversity_ov(ps,'Phylum')
all.div.plots_gen <-alpha_diversity_ov(ps,'Genus')
pop_full_barplot <- barplot_ov(ps,'Phylum',full=T,abd=0.3)
pop_gen_barplot <- barplot_ov(ps,'Genus',full=T,abd=0.3)
beta_pop_plot <- beta_diversity_pop(ps,clade="None",loc="",prev_th=0.1,type_dist='bray')
beta_phy_pop_plot <- beta_diversity_pop(ps,clade="Phylum",loc="",prev_th=0.1,type_dist='bray')
beta_gen_pop_plot <- beta_diversity_pop(ps,clade="Genus",loc="",prev_th=0.1,type_dist='bray')
beta_ov_plot <- beta_diversity_ov(ps,clade="None",loc="",prev_th=0.1,type_dist='bray')
beta_ov_phy_plot <- beta_diversity_ov(ps,clade="Phylum",loc="",prev_th=0.1,type_dist='bray')
beta_ov_gen_plot <- beta_diversity_ov(ps,clade="Genus",loc="",prev_th=0.1,type_dist='bray')

ps2.venn <- merge_samples(filtered_ps, 'Population', fun = sum)
venn_obj <- as.data.frame(t(otu_table(ps2.venn)))
venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, 1, 0),
                          USE.NAMES = T)
rownames(venn_obj.binary) <- rownames(venn_obj)
venn_obj.binary <- as.data.frame(venn_obj.binary)
upset_order <- colnames(venn_obj.binary)
shared_ASV_plot_pop <- upset(venn_obj.binary, nsets = 6,
                         sets = rev(upset_order),
                         mainbar.y.label = 'Shared Genera',
                         sets.x.label = 'Genera per Group',
                         keep.order = T,
                         order.by = 'freq', 
                         sets.bar.color = rev(primary_color_list),
                         text.scale = 1.6)
png('./pop_genus_shared_genus.png', width=6.5,units="in", height=5, res=300)
print(shared_ASV_plot_pop)
dev.off()

ps2.venn <- merge_samples(filtered_ps, 'Overwintered', fun = sum)
venn_obj <- as.data.frame(t(otu_table(ps2.venn)))
venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, 1, 0),
                          USE.NAMES = T)
rownames(venn_obj.binary) <- rownames(venn_obj)
venn_obj.binary <- as.data.frame(venn_obj.binary)
upset_order <- colnames(venn_obj.binary)
shared_ASV_plot_ov <- upset(venn_obj.binary, nsets = 6,
                         sets = rev(upset_order),
                         mainbar.y.label = 'Shared Genera',
                         sets.x.label = 'Genera per Group',
                         keep.order = T,
                         order.by = 'freq', 
                         sets.bar.color = c("#1B9E77","#D95F02"),
                         text.scale = 1.6)
png('./ov_genus_shared_genus.png', width=6.5,units="in", height=5, res=300)
print(shared_ASV_plot_ov)
dev.off()

####### overwintered population by population analyses ######

dev.off()
ps.TN <- prune_samples(sample_data(ps)$Population == "TN", ps)
ps.TN <- prune_taxa(taxa_sums(ps.TN) > 0, ps.TN)
setwd(paste('./','TN',sep=''))
barplot_ov(ps.TN,'Phylum','TN')
barplot_ov(ps.TN,'Genus','TN')
TN.div.plots <-alpha_diversity_ov(ps.TN,'Phylum','TN')
shan.TN.ov.phy <- TN.div.plots[[1]]
chao.TN.ov.phy <- TN.div.plots[[2]]
TN.div.plots <- alpha_diversity_ov(ps.TN,'Genus','TN')
shan.TN.ov.gen <- TN.div.plots[[1]]
chao.TN.ov.gen <- TN.div.plots[[2]]
TN.sigtab <- DESeq_ov.pop(ps,'TN',clade="Genus")
TN.deseq.plot <- TN.sigtab[[2]]
TN.sigtab <- TN.sigtab[[1]]
TN.sig.asv <- rownames(TN.sigtab)
TN.sig.proportion <- dim(TN.sigtab)[1]/dim(otu_table(ps.TN))[2]
ASVtab.TN.sig <- otu_table(ps)[,which(colnames(otu_table(ps.TN)) %in% TN.sig.asv)]
OTU1 <- as(ASVtab.TN.sig, "matrix")
OTUdf <- data.frame(t(OTU1)) # taxa rows, samples columns
OTUbiom <- make_biom(OTUdf)
write_biom(OTUbiom, 'OTU_table.biom')
TN.asv.seq <- refseq.df[which(rownames(refseq.df) %in% TN.sig.asv),]
TN.dna <- Biostrings::DNAStringSet(TN.asv.seq)
names(TN.dna) <- TN.sig.asv
TN.dna%>%
  Biostrings::writeXStringSet("./asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")
pantoea_check(ps.TN)


ps.BS <- prune_samples(sample_data(ps)$Population == "BS", ps)
ps.BS <- prune_taxa(taxa_sums(ps.BS) > 0, ps.BS)
setwd(paste('../','BS',sep=''))
barplot_ov(ps.BS,'Phylum','BS')
barplot_ov(ps.BS,'Genus','BS')
BS.div.plots <-alpha_diversity_ov(ps.BS,'Phylum','BS')
shan.BS.ov.phy <- BS.div.plots[[1]]
chao.BS.ov.phy <- BS.div.plots[[2]]
BS.div.plots <- alpha_diversity_ov(ps.BS,'Genus','BS')
shan.BS.ov.gen <- BS.div.plots[[1]]
chao.BS.ov.gen <- BS.div.plots[[2]]
BS.sigtab <- DESeq_ov.pop(ps,'BS',clade="Genus")
BS.deseq.plot <- BS.sigtab[[2]]
BS.sigtab <- BS.sigtab[[1]]
BS.sig.asv <- rownames(BS.sigtab)
BS.sig.proportion <- dim(BS.sigtab)[1]/dim(otu_table(ps.BS))[2]
ASVtab.BS.sig <- otu_table(ps)[,which(colnames(otu_table(ps.BS)) %in% BS.sig.asv)]
OTU1 <- as(ASVtab.BS.sig, "matrix")
OTUdf <- data.frame(t(OTU1)) # taxa rows, samples columns
OTUbiom <- make_biom(OTUdf)
write_biom(OTUbiom, 'OTU_table.biom')
BS.asv.seq <- refseq.df[which(rownames(refseq.df) %in% BS.sig.asv),]
BS.dna <- Biostrings::DNAStringSet(BS.asv.seq)
names(BS.dna) <- BS.sig.asv
BS.dna%>%
  Biostrings::writeXStringSet("./asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")
pantoea_check(ps.BS)

ps.DN <- prune_samples(sample_data(ps)$Population == 'DN', ps)
ps.DN <- prune_taxa(taxa_sums(ps.DN) > 0, ps.DN)
setwd(paste('../','DN',sep=''))
barplot_ov(ps.DN,'Phylum','DN')
barplot_ov(ps.DN,'Genus','DN')
DN.div.plots <-alpha_diversity_ov(ps.DN,'Phylum','DN')
shan.DN.ov.phy <- DN.div.plots[[1]]
chao.DN.ov.phy <- DN.div.plots[[2]]
DN.div.plots <- alpha_diversity_ov(ps.DN,'Genus','DN')
shan.DN.ov.gen <- DN.div.plots[[1]]
chao.DN.ov.gen <- DN.div.plots[[2]]
DN.sigtab <- DESeq_ov.pop(ps,'DN',clade="Genus")
DN.deseq.plot <- DN.sigtab[[2]]
DN.sigtab <- DN.sigtab[[1]]
DN.sig.asv <- rownames(DN.sigtab)
DN.sig.proportion <- dim(DN.sigtab)[1]/dim(otu_table(ps.DN))[2]
ASVtab.DN.sig <- otu_table(ps)[,which(colnames(otu_table(ps.DN)) %in% DN.sig.asv)]
OTU1 <- as(ASVtab.DN.sig, "matrix")
OTUdf <- data.frame(t(OTU1)) # taxa rows, samples columns
OTUbiom <- make_biom(OTUdf)
write_biom(OTUbiom, 'OTU_table.biom')
DN.asv.seq <- refseq.df[which(rownames(refseq.df) %in% DN.sig.asv),]
DN.dna <- Biostrings::DNAStringSet(DN.asv.seq)
names(DN.dna) <- DN.sig.asv
DN.dna%>%
  Biostrings::writeXStringSet("./asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")
pantoea_check(ps.DN)

ps.SM <- prune_samples(sample_data(ps)$Population == 'SM', ps)
ps.SM <- prune_taxa(taxa_sums(ps.SM) > 0, ps.SM)
setwd(paste('../','SM',sep=''))
barplot_ov(ps.SM,'Phylum','SM')
barplot_ov(ps.SM,'Genus','SM')
SM.div.plots <-alpha_diversity_ov(ps.SM,'Phylum','SM')
shan.SM.ov.phy <- SM.div.plots[[1]]
chao.SM.ov.phy <- SM.div.plots[[2]]
SM.div.plots <- alpha_diversity_ov(ps.SM,'Genus','SM')
shan.SM.ov.gen <- SM.div.plots[[1]]
chao.SM.ov.gen <- SM.div.plots[[2]]
SM.sigtab <- DESeq_ov.pop(ps,'SM',clade="Genus")
SM.deseq.plot <- SM.sigtab[[2]]
SM.sigtab <- SM.sigtab[[1]]
SM.sig.asv <- rownames(SM.sigtab)
SM.sig.proportion <- dim(SM.sigtab)[1]/dim(otu_table(ps.SM))[2]
ASVtab.SM.sig <- otu_table(ps)[,which(colnames(otu_table(ps.SM)) %in% SM.sig.asv)]
OTU1 <- as(ASVtab.SM.sig, "matrix")
OTUdf <- data.frame(t(OTU1)) # taxa rows, samples columns
OTUbiom <- make_biom(OTUdf)
write_biom(OTUbiom, 'OTU_table.biom')
SM.asv.seq <- refseq.df[which(rownames(refseq.df) %in% SM.sig.asv),]
SM.dna <- Biostrings::DNAStringSet(SM.asv.seq)
names(SM.dna) <- SM.sig.asv
SM.dna%>%
  Biostrings::writeXStringSet("./asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")
pantoea_check(ps.SM)


setwd('..')

ps <- tax_glom(ps,"Genus",NArm=F)
ps <- prune_samples(sample_data(ps)$Population != 'VIG', ps)
# Prevalence filtering to 10% prevalence
prevThreshold <- nsamples(ps) * 0.10

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevThreshold)]
removeTaxa <- rownames(prevdf)[(prevdf$Prevalence < prevThreshold)]
ps2 <- prune_taxa(keepTaxa, ps)
ps2.comp <- transform_sample_counts(ps2,  
                                    function(x) x/sum(x) * 100)
# Get matching metadata & ASV table

ord.meta <- data.frame(sample_data(ps2))
ord.asvs <- otu_table(ps2.comp)

# Generate nmds & pcoa ordinations (Bray-Curtis dissimilarity)
bc.nmds <- metaMDS(ord.asvs, 
                   autotransform = F,
                   trymax = 100)
bc.dist <- metaMDSredist(bc.nmds)
bc.envfit.nmds <- envfit(bc.nmds ~ Population,
                         data = ord.meta,
                         display = 'sites',
                         na.rm = T)

# Extract fit data
# NMDS data
# Use ggvegan to extract the nmds data for plotting as ggplot object
bc.envfit.nmds.fort <- fortify(bc.envfit.nmds)
centroids.nmds <- subset(bc.envfit.nmds.fort, Type == 'Centroid')
centroids.nmds$Label <- c('SM', 'TN', "BS", "DN")
vectors.nmds <- subset(bc.envfit.nmds.fort, Type == 'Vector')

# Get ellipse data
# Ellipse function
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

bc.nmds.df <- data.frame(MDS1 = bc.nmds$points[,1],
                         MDS2 = bc.nmds$points[,2],
                         Population = ord.meta$Population)

plot.new()
bc.nmds.ellipse <- ordiellipse(bc.nmds, 
                               ord.meta$Population, 
                               kind = "sd", 
                               conf = 0.95,
                               object = T)
bc.nmds.ellipse.df <- data.frame()
ord.meta$Population <- as.factor(ord.meta$Population)
for(tp in levels(ord.meta$Population)){
  bc.nmds.ellipse.df <- rbind(bc.nmds.ellipse.df, 
                              cbind(as.data.frame(with(bc.nmds.df[bc.nmds.df$Population==tp,],
                                                       veganCovEllipse(bc.nmds.ellipse[[tp]]$cov,
                                                                       bc.nmds.ellipse[[tp]]$center,
                                                                       bc.nmds.ellipse[[tp]]$scale)))
                                    ,Population=tp))
}

# vectors.pcoa <- subset(bc.envfit.pcoa.fort, Type == 'Vector')

ord.bray.nmds.plot <- plot_ordination(
  ps2.comp,
  bc.nmds,
  color = 'Population') +
  theme_cowplot() + 
  background_grid(major = 'xy', minor = 'none') +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = primary_color_list)

ord.bray.nmds.plot +
  geom_path(data = bc.nmds.ellipse.df,
            aes(x = NMDS1, y = NMDS2, color = Population),
            size = 1, 
            linetype = 2,
            inherit.aes = F) 

bc.adonis <- adonis2(bc.dist ~ Population,
                     data = ord.meta)
bc.adonis

bc.disper <- betadisper(bc.dist, ord.meta$Population,type="centroid")
TukeyHSD(bc.disper)
permutest(bc.disper)
pairwise.adonis(
  otu_table(ps),
  sample_data(ps)$Population,
  sim.method = "bray",
  p.adjust.m = "none"
)
png('betadisp_ov_tot.png',width=(3.25*1200/72), height=(3.25*1200/72), res=1200)
plot(bc.disper, hull=FALSE, ellipse=TRUE)
dev.off()

ps1 <- prune_taxa(taxa_sums(ps2) > 0, ps2)
save.image(file = "lefse_preparation.RData")

######## beta div pca population #############
# based on https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

#centred log ratio transform
(ps_clr <- microbiome::transform(ps.gen, "clr"))    
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps.gen, "NMDS")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps.gen, ord_clr, type="samples", color="Population") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Population), linetype = 2)
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps.gen, method = "bray") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps.gen)$Population)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps.gen)$Population)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
df.dispr <- as.data.frame(cbind(dispr$distances,sample_data(ps.gen)$Population))
colnames(df.dispr)<-c('distances','Population')
df.dispr$distances <- as.numeric(df.dispr$distances)
wilcox.centr <- df.dispr %>%
  wilcox_test(distances ~ Population) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")
# df.dispr$Population <- factor(df.dispr$Population, levels = c("BS","TN","SM","DN"))
Bdist_distance_pop <- ggboxplot(df.dispr, x = "Population", y = "distances",
         color = "black", palette = primary_color_list, fill = "Population") +
  stat_pvalue_manual(wilcox.centr, label = "p.adj.signif", y.position = (max(df.dispr$distances))+0.1, step.increase = 0.1) +
   ylab('Distance from centroids') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")
permutest(dispr)
#Generate distances
ord_unifrac <- ordinate(ps, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps, ord_unifrac, color = "Population") + geom_point(size = 2)
b <- plot_ordination(ps, ord_unifrac_un, color = "Population") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

######## beta div pca population #############
# based on https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

#centred log ratio transform
ps.phy <- tax_glom(ps,"Phylum",NArm=F)
(ps_clr <- microbiome::transform(ps.phy, "clr"))    
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps.phy, "NMDS")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps.phy, ord_clr, type="samples", color="Population") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Population), linetype = 2)
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps.phy, method = "bray") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps.phy)$Population)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps.phy)$Population)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
df.dispr <- as.data.frame(cbind(dispr$distances,sample_data(ps.phy)$Population))
colnames(df.dispr)<-c('distances','Population')
df.dispr$distances <- as.numeric(df.dispr$distances)
wilcox.centr <- df.dispr %>%
  wilcox_test(distances ~ Population) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")
# df.dispr$Population <- factor(df.dispr$Population, levels = c("BS","TN","SM","DN"))
Bdist_distance_pop_phy <- ggboxplot(df.dispr, x = "Population", y = "distances",
                                color = "black", palette = primary_color_list, fill = "Population") +
  stat_pvalue_manual(wilcox.centr, label = "p.adj.signif", y.position = (max(df.dispr$distances))+0.1, step.increase = 0.1) +
  ylab('Distance from centroids') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")
permutest(dispr)
#Generate distances
ord_unifrac <- ordinate(ps, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps, ord_unifrac, color = "Population") + geom_point(size = 2)
b <- plot_ordination(ps, ord_unifrac_un, color = "Population") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

# based on https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#centred log ratio transform
ps.prop<- transform_sample_counts(ps, function(x) x/sum(x)*100)
(ps_clr <- microbiome::transform(ps.prop, "clr"))    
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
data$Treatment <- factor(data$Treatment, levels=unique(data$Treatment))
phyloseq::plot_ordination(ps.prop, ord_clr, type="samples", color="Population") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Population), linetype = 2)+
  scale_color_manual(values=primary_color_list,breaks = c('BS','TN','SM','DN'))
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "bray") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Population)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Population)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr)
#Generate distances
ord_unifrac <- ordinate(ps.prop, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps.prop, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps.prop, ord_unifrac, color = "Population") + geom_point(size = 2)
b <- plot_ordination(ps.prop, ord_unifrac_un, color = "Population") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

######## beta div pca #############
# based on https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#centred log ratio transform
(ps_clr <- microbiome::transform(ps.phy, "clr"))    
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps.phy, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps.phy, ord_clr, type="samples", color="Overwintered") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Overwintered), linetype = 2)
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps.phy, method = "bray") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps.phy)$Overwintered)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps.phy)$Overwintered)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
df.dispr <- as.data.frame(cbind(dispr$distances,sample_data(ps.phy)$Overwintered))
colnames(df.dispr)<-c('distances','Overwintered')
df.dispr$distances <- as.numeric(df.dispr$distances)
wilcox.centr <- df.dispr %>%
  wilcox_test(distances ~ Overwintered) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")
df.dispr$Overwintered <- factor(df.dispr$Overwintered)
Bdist_distance_ov_phy <- ggboxplot(df.dispr, x = "Overwintered", y = "distances",
                                color = "black", palette = primary_color_list, fill = "Overwintered",breaks = c('No', 'Yes')) +
  stat_pvalue_manual(wilcox.centr, label = "p.adj.signif", y.position = (max(df.dispr$distances))+0.1, step.increase = 0.1) +
  ylab('Distance from centroids') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")
permutest(dispr)
#Generate distances
ord_unifrac <- ordinate(ps, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps, ord_unifrac, color = "Overwintered") + geom_point(size = 2)
b <- plot_ordination(ps, ord_unifrac_un, color = "Overwintered") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

######## beta div pca genus #############
# based on https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#centred log ratio transform
ps_gen <- ps.gen
ps_gen_prop <- transform_sample_counts(ps_gen, function(x) x/sum(x)*100)
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_gen_prop, method = "bray") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_gen)$Overwintered)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_gen)$Overwintered)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
df.dispr <- as.data.frame(cbind(dispr$distances,sample_data(ps_gen)$Overwintered))
colnames(df.dispr)<-c('distances','Overwintered')
df.dispr$distances <- as.numeric(df.dispr$distances)
wilcox.centr <- df.dispr %>%
  wilcox_test(distances ~ Overwintered) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")
df.dispr$Overwintered <- factor(df.dispr$Overwintered)
Bdist_gen_distance_ov <- ggboxplot(df.dispr, x = "Overwintered", y = "distances",
                               color = "black", palette = primary_color_list, fill = "Overwintered",breaks = c('No', 'Yes')) +
  stat_pvalue_manual(wilcox.centr, label = "p.adj.signif", y.position = (max(df.dispr$distances))+0.1, step.increase = 0.1) +
  ylab('Distance from centroids') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")
permutest(dispr)

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_gen)$Population)
#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_gen)$Population)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
df.dispr <- as.data.frame(cbind(dispr$distances,sample_data(ps_gen)$Population))
colnames(df.dispr)<-c('distances','Population')
df.dispr$distances <- as.numeric(df.dispr$distances)
wilcox.centr <- df.dispr %>%
  wilcox_test(distances ~ Population) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")
df.dispr$Population <- factor(df.dispr$Population, levels=c('BS', 'TN', "SM", 'DN'))
Bdist_gen_distance_pop <- ggboxplot(df.dispr, x = "Population", y = "distances",
                                   color = "black", palette = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), fill = "Population",breaks = c('SM', 'TN', "BS", 'DN')) +
  stat_pvalue_manual(wilcox.centr, label = "p.adj.signif", y.position = (max(df.dispr$distances))+0.1, step.increase = 0.1) +
  ylab('Distance from centroids') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")
permutest(dispr)

####### Nosema maddoxxii part #######
sum_nosema <- c(1/16,6/15,0,0,0,1/15,0,0)
pop_nosema <- c('TN','BS','SM','DN','TN','BS','SM','DN')
ov_nosema <- c(rep('No',4),rep('Yes',4))
df_nosema <- as.data.frame(cbind(pop_nosema,as.numeric(sum_nosema),ov_nosema))
colnames(df_nosema) <- c('Population','% Positive','Overwintered')
df_nosema$Population <- as.factor(df_nosema$Population)
df_nosema$Overwintered <- as.factor(df_nosema$Overwintered)
df_nosema$`% Positive` <- as.numeric(df_nosema$`% Positive`)

bp = ggplot(data=df_nosema, aes(x=Population, y=`% Positive`,fill=Overwintered)) +
  geom_bar(stat="identity",position = position_stack(reverse = TRUE))+
  scale_fill_manual(values=primary_color_list,breaks = c('Yes','No'))

png('./nosema_barplot.png', width=6.5,units="in", height=5, res=1200)
print(bp)
dev.off()

######## paper plots section #######

load("G:/Il mio Drive/PhD.AES.Dsuz.Riccardo.Piccinno/H_halys/metaBUG/Sequenze_MiSeq/complete_analysis_figures_paper.RData")
library(png)
shared_ASV_plot_pop <- rasterGrob(readPNG("pop_genus_shared_genus.png"))
shared_ASV_plot_ov <- rasterGrob(readPNG("ov_genus_shared_genus.png"))

# Figure 3 - population plots
top_row <- plot_grid(pop_gen_barplot+guides(x='none')+ theme(text = element_text(size = 12)),
                     shared_ASV_plot_pop,
                     labels = c('A', 'B'), label_size = 12)
mid_row <- plot_grid(beta_gen_pop_plot + geom_point(size=0.01) + scale_shape_manual(values = c(1,2,3,4), breaks=c('BS','TN','SM','DN')) + theme(text = element_text(size = 16)), 
                     Bdist_distance_pop +theme(legend.position = 'none',text = element_text(size = 12)),
                     labels = c('A', 'B'), label_size = 12, rel_widths = c(2,1.7), ncol=2, nrow=1)
bot_row <- plot_grid(pop.div.plots_gen[[1]] + theme(text = element_text(size = 12)),
                     pop.div.plots_gen[[2]] +theme(legend.position = 'left',text = element_text(size = 12)), 
                     labels = c('C',' D'), label_size = 12, rel_widths = c(1.7,2), ncol=2, nrow=1)
png('Figure3.png', width=10,units="in", height=13, res=300)
plot_grid(mid_row,bot_row,nrow=2, ncol=1)
dev.off()

# Figure 4 - overwintering plots genus
top_row <- plot_grid(beta_ov_gen_plot + scale_shape_manual(values = c(1,1,1,1), breaks=c('No','Yes')) + theme(text = element_text(size = 16)),
                     Bdist_gen_distance_ov  + theme(text = element_text(size = 16)),
                     labels = c('A', 'B'), label_size = 12)
bot_row <- plot_grid(all.div.plots_gen[[1]] + theme(legend.position = 'right',text = element_text(size = 16)),
                     all.div.plots_gen[[2]] +theme(text = element_text(size = 16)), 
                     labels = c('C', 'D'), label_size = 12, rel_widths = c(1.5,1.5), ncol=2, nrow=1)
png('Figure4.png', width=13,units="in", height=10, res=300)
plot_grid(top_row,bot_row,nrow=2, ncol=1)
dev.off()

# Figure 5 - overwinter pop wise + deseq
png('Figure5.png', width=10,units="in", height=13, res=300)
plot_grid(shan.TN.ov.gen + theme(text = element_text(size = 16)),
          chao.TN.ov.gen + theme(text = element_text(size = 16),legend.position = 'none'),
          TN.deseq.plot$gtable,
          shan.BS.ov.gen + theme(text = element_text(size = 16)),
          chao.BS.ov.gen + theme(text = element_text(size = 16),legend.position = 'none'),
          BS.deseq.plot$gtable,
          shan.SM.ov.gen + theme(text = element_text(size = 16)),
          chao.SM.ov.gen + theme(text = element_text(size = 16),legend.position = 'none'),
          SM.deseq.plot$gtable,
          shan.DN.ov.gen + theme(text = element_text(size = 16)),
          chao.DN.ov.gen + theme(text = element_text(size = 16),legend.position = 'none'),
          DN.deseq.plot$gtable,
          labels = c('A', 'B', 'C', 'D' , 'E', 'F', 'G', 'H', "I", "L", "M", "N"), label_size = 10,
          nrow=4,ncol=3,rel_widths = c(0.75,0.75,2,0.75,0.75,2,0.75,0.75,2,0.75,0.75,2),
          scale = c(.94,.94,.94,.94,.94,.94,.94,.94,.94,.94,.94,.94))
dev.off()

# Figure S1 - population plots
top_row <- plot_grid(pop_full_barplot+guides(x='none')+ theme(text = element_text(size = 16)), 
                     beta_phy_pop_plot + geom_point(size=0.01) + scale_shape_manual(values = c(1,2,3,4), breaks=c('BS','TN','SM','DN')) + theme(text = element_text(size = 16)), 
                     labels = c('A', 'B'), label_size = 12)
bot_row <- plot_grid(pop.div.plots_phy[[1]] + theme(text = element_text(size = 16)),
                     pop.div.plots_phy[[2]] +theme(legend.position = 'none',text = element_text(size = 16)), 
                     Bdist_distance_pop_phy + theme(text = element_text(size = 16)),
                     labels = c('C', 'D',' E'), label_size = 12, rel_widths = c(1.5,1.5,2), ncol=3, nrow=1)
png('FigureS1.png', width=13,units="in", height=10, res=300)
plot_grid(top_row,bot_row,nrow=2, ncol=1)
dev.off()

# Figure S2 - overwintering plots phylum
top_row <- plot_grid(beta_ov_phy_plot + scale_shape_manual(values = c(1,2), breaks=c('No','Yes')) + theme(text = element_text(size = 16)),
                     Bdist_distance_ov_phy + theme(text = element_text(size = 16)),
                     labels = c('A', 'B'), label_size = 12)
bot_row <- plot_grid(all.div.plots_phy[[1]] + theme(legend.position = 'right',text = element_text(size = 16)),
                     all.div.plots_phy[[2]] +theme(text = element_text(size = 16)), 
                     labels = c('C', 'D'), label_size = 12, rel_widths = c(1.5,1.5), ncol=2, nrow=1)
png('FigureS2.png', width=13,units="in", height=10, res=300)
plot_grid(top_row,bot_row,nrow=2, ncol=1)
dev.off()

# Figure S3 - overwinter pop wise
png('FigureS3.png', width=13,units="in", height=10, res=300)
plot_grid(shan.TN.ov.phy + theme(text = element_text(size = 16)),
          chao.TN.ov.phy + theme(text = element_text(size = 16),legend.position = 'none'),
          shan.BS.ov.phy + theme(text = element_text(size = 16)),
          chao.BS.ov.phy + theme(text = element_text(size = 16)),
          shan.SM.ov.phy + theme(text = element_text(size = 16)),
          chao.SM.ov.phy + theme(text = element_text(size = 16),legend.position = 'none'),
          shan.DN.ov.phy + theme(text = element_text(size = 16)),
          chao.DN.ov.phy + theme(text = element_text(size = 16)),
          labels = c('A', 'B', 'C', 'D' , 'E', 'F', 'G', 'H'), label_size = 12,
          nrow=2,ncol=4,rel_widths = c(1.25,1.25,1.25,1.75,1.25,1.25,1.25,1.75))
dev.off()






