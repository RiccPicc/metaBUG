barplot_ov <- function(ps,clade="None",loc="",full=F,abd=0.15){
  if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
  ps.prop <- transform_sample_counts(ps, function(x) x/sum(x)*100)
  ps2.bar <- ps.prop
  
  my_plot_bar <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                           facet_grid = NULL) {
    mdf = psmelt(physeq)
    p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
    p = p + geom_bar(stat = "identity", position = "stack") + xlab('') + labs(x='') + scale_fill_brewer(palette="Set2")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 1))
    if (!is.null(facet_grid)) {
      p <- p + facet_grid(facet_grid)
    }
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    return(p)
  }
  
  # Merged phyloseq plot (100% total abundance)
  ps2.bar.merge <- merge_samples(ps2.bar, group = 'Overwintered')
  ## Warning in asMethod(object): NAs introduced by coercion
  sample_data(ps2.bar.merge)$Overwintered <- sample_names(ps2.bar.merge)
  ps2.bar.merge.comp <- transform_sample_counts(ps2.bar.merge,  
                                                function(x) x/sum(x) * 100)
  ps2.bp <- tax_glom(ps2.bar.merge.comp, clade, NArm = F)
  cabd <- sort(taxa_sums(ps2.bp)/nsamples(ps2.bp),
               decreasing = T)
  # Lets remove some low abundance phyla
  renamed.clade <- names(which(cabd < abd))
  tax_table(ps2.bp)[renamed.clade, clade] <- paste('< ',abd,'% Abd',sep="")
  ps2.bp <- tax_glom(ps2.bp, clade, NArm = F)
  cabd <- sort(taxa_sums(ps2.bp)/nsamples(ps2.bp),
               decreasing = T)
  taxa_count <- length(unique(tax_table(ps2.bp)[,clade]))
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  
  taxa_order <- as.character(tax_table(ps2.bp)[names(rev(cabd)), clade])
  taxa_order <- unique(taxa_order)
  sample_order <- sample_names(ps2.bp)
  bp <- plot_bar(ps2.bp, fill = clade) +
    geom_bar(stat='identity',
             position='stack') +
    scale_fill_manual(values = rev(getPalette(taxa_count))) + 
    scale_color_manual(name = clade,
                       values = rev(getPalette(taxa_count))) +
    ylim(c(0,101)) +
    theme_cowplot() +
    xlab('Overwintered')
  bp$data[,clade] <- factor(bp$data[,clade], 
                                      levels = taxa_order)
  bp$data[,'Sample'] <- factor(bp$data[,'Sample'], 
                                      levels = sample_order)
  
  png(paste('./ov_',clade,'_',loc,'_barplot.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(bp)
  dev.off()
  
  ps.prop.renamed <- ps.prop
  cabd <- sort(taxa_sums(ps.prop.renamed)/nsamples(ps.prop.renamed),
                                  decreasing = T)
  # Lets remove some low abundance phyla
  renamed.clade <- names(which(cabd < abd))
  tax_table(ps.prop.renamed)[renamed.clade, clade] <- paste('< ',abd,'% Abd',sep="")
  print(abd)
  png(paste('./pop_phylobars_full_',clade,'.png',sep=''), width=9.5,units="in", height=5, res=1200)
  print(my_plot_bar(ps.prop.renamed, fill=clade, x="Subject") + facet_wrap(~Population, scales="free_x", nrow=1)+guides(x='none')+ theme(text = element_text(size = 16)))
  dev.off()
  
  png(paste('./ov_phylobars_full_',clade,'.png',sep=''), width=9.5,units="in", height=5, res=1200)
  print(my_plot_bar(ps.prop.renamed, fill=clade, x="Subject") + facet_wrap(~Overwintered, scales="free_x", nrow=1)+guides(x='none')+ theme(text = element_text(size = 16)))
  dev.off()
  
  # if (full==T){
  #   # ampvis2
  #   #Combine OTU abundance table and taxonomy table from the phyloseq object "my_phyloseq_object":
  #   av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(ps2.bar)@.Data)),
  #                              t(phyloseq::otu_table(ps2.bar)@.Data),
  #                              phyloseq::tax_table(ps2.bar)@.Data,
  #                              check.names = F
  #   )
  #   
  #   #Extract metadata from the phyloseq object:
  #   av2_metadata <- data.frame(phyloseq::sample_data(ps2.bar), 
  #                              check.names = F
  #   )
  #   
  #   av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)
  #   
  #   #Load the data with amp_load:
  #   av2_obj <- amp_load(av2_otutable, av2_metadata)
  #   av2_box <- amp_boxplot(av2_obj, group_by = 'Overwintered',
  #                          sort_by = 'mean', tax_show = 10,
  #                          tax_aggregate = 'Phylum', 
  #                          normalise = F) +
  #     scale_color_manual(values = c("#D95F02","#1B9E77"))
  #   png(paste('./ov_',clade,'_',loc,'_barplot_full_phylum.png',sep=''), width=6.5,units="in", height=5, res=1200)
  #   print(av2_box)
  #   dev.off()
  #   
  #   av2_box <- amp_boxplot(av2_obj, group_by = 'Overwintered',
  #                          sort_by = 'mean', tax_show = 10,
  #                          tax_aggregate = 'Genus',
  #                          normalise = F) +
  #     scale_color_manual(values = c("#D95F02","#1B9E77"))
  #   png(paste('./ov_',clade,'_',loc,'_barplot_full_genus.png',sep=''), width=6.5,units="in", height=5, res=1200)
  #   print(av2_box)
  #   dev.off()
  #   
  #   #Load the data with amp_load:
  #   av2_obj <- amp_load(av2_otutable, av2_metadata)
  #   av2_box <- amp_boxplot(av2_obj, group_by = 'Population',
  #                          sort_by = 'mean', tax_show = 10,
  #                          tax_aggregate = 'Phylum', 
  #                          normalise = F) +
  #     scale_color_manual(values = rev(primary_color_list))
  #   png(paste('./pop_',clade,'_',loc,'_barplot_full_phylum.png',sep=''), width=6.5,units="in", height=5, res=1200)
  #   print(av2_box)
  #   dev.off()
  #   
  #   av2_box <- amp_boxplot(av2_obj, group_by = 'Population',
  #                          sort_by = 'mean', tax_show = 10,
  #                          tax_aggregate = 'Genus',
  #                          normalise = F) +
  #     scale_color_manual(values = rev(primary_color_list))
  #   png(paste('./pop_',clade,'_',loc,'_barplot_full_genus.png',sep=''), width=6.5,units="in", height=5, res=1200)
  #   print(av2_box)
  #   dev.off()
  #   
  #   # Use UpSetR to generate upset plots (Venn diagram alternative)
  #   ps2.venn <- merge_samples(ps2.bar, 'Overwintered', fun = sum)
  #   ## Warning in asMethod(object): NAs introduced by coercion
  #   venn_obj <- as.data.frame(t(otu_table(ps2.venn)))
  #   venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, 1, 0),
  #                             USE.NAMES = T)
  #   rownames(venn_obj.binary) <- rownames(venn_obj)
  #   venn_obj.binary <- as.data.frame(venn_obj.binary)
  #   upset_order <- colnames(venn_obj.binary)
  #   shared_ASV_plot <- upset(venn_obj.binary, nsets = 6,
  #                            sets = rev(upset_order),
  #                            mainbar.y.label = 'Shared ASVs',
  #                            sets.x.label = 'ASVs per Group',
  #                            keep.order = T,
  #                            order.by = 'freq', sets.bar.color = c("#1B9E77","#D95F02"))
  #   png(paste('./ov_',clade,'_',loc,'_shared_asv.png',sep=''), width=6.5,units="in", height=5, res=1200)
  #   print(shared_ASV_plot)
  #   dev.off()
  #   
  #   ps2.venn <- merge_samples(ps, 'Population', fun = sum)
  #   ## Warning in asMethod(object): NAs introduced by coercion
  #   venn_obj <- as.data.frame(t(otu_table(ps2.venn)))
  #   venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, 1, 0),
  #                             USE.NAMES = T)
  #   rownames(venn_obj.binary) <- rownames(venn_obj)
  #   venn_obj.binary <- as.data.frame(venn_obj.binary)
  #   upset_order <- colnames(venn_obj.binary)
  #   shared_ASV_plot <- upset(venn_obj.binary, nsets = 6,
  #                            sets = rev(upset_order),
  #                            mainbar.y.label = 'Shared ASVs',
  #                            sets.x.label = 'ASVs per Group',
  #                            keep.order = T,
  #                            order.by = 'freq', sets.bar.color = rev(primary_color_list))
  #   png(paste('./pop_',clade,'_',loc,'_shared_asv.png',sep=''), width=6.5,units="in", height=5, res=1200)
  #   print(shared_ASV_plot)
  #   dev.off()
    
    return(my_plot_bar(ps.prop.renamed, fill=clade, x="Subject") + facet_wrap(~Population, scales="free_x", nrow=1))
  }
alpha_diversity_pop <- function(ps,clade="None",loc="") {
if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
richness.rare <- cbind(estimate_richness(ps, 
                                         measures = c('Shannon', 'Chao1')),
                       sample_data(ps)$Population)
colnames(richness.rare) <- c('Chao1','se.chao1', 'Shannon', 'Population')
richness.rare$Labels <- rownames(richness.rare)
ps.prop <- transform_sample_counts(ps, function(x) x/sum(x)*100)

ad.test.df <- richness.rare[,c('Shannon', 'Chao1')]
ad.test.df <- cbind(ad.test.df,
                    sample_data(ps.prop))

ad.wilcox.shannon <- ad.test.df %>%
  wilcox_test(Shannon ~ Population) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")

ad.wilcon.Chao1 <- ad.test.df %>%
  wilcox_test(Chao1 ~ Population) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.adj.signif != "ns")

write_tsv(ad.test.df %>%
        wilcox_test(Shannon ~ Population) %>%
        adjust_pvalue() %>%
        add_significance(), paste(clade, "Shannon_test.tsv"))
write_tsv(ad.test.df %>%
        wilcox_test(Chao1 ~ Population) %>%
        adjust_pvalue() %>%
        add_significance(), paste(clade,"Chao1_test.tsv"))

shannon.plot.sig <- ggboxplot(ad.test.df, x = "Population", y = "Shannon",
                              color = "black", palette = primary_color_list, fill = "Population") +
  stat_pvalue_manual(ad.wilcox.shannon, label = "p.adj.signif", y.position = (max(ad.test.df$Shannon))+0.1, step.increase = 0.1) + 
  ylab('Shannon Index') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")

Chao1.plot.sig <- ggboxplot(ad.test.df, x = "Population", y = "Chao1",
                            color = "black", palette = primary_color_list, fill = "Population") + 
  stat_pvalue_manual(ad.wilcon.Chao1, label = "p.adj.signif", y.position = (max(ad.test.df$Chao1))+0.1, step.increase = 0.1) + 
  ylab('Chao1 Index') +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) + 
  grids(linetype = "dashed")


diversity.plots.sig <- plot_grid(shannon.plot.sig, Chao1.plot.sig, 
                                 labels = c('A', 'B'), align = 'h',rel_widths = c(1.5, 2))  
png(paste('./pop_',clade,'_',loc,'_alpha_div.png',sep=''), width=6.5,units="in", height=5, res=1200)
print(diversity.plots.sig)
dev.off()
print(ad.wilcox.shannon)

return(list(shannon.plot.sig,Chao1.plot.sig))
}
alpha_diversity_ov <- function(ps,clade="None",loc=""){
  if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
  richness.rare <- cbind(estimate_richness(ps, 
                                           measures = c('Shannon', 'Chao1')),
                         sample_data(ps)$Overwintered)
  colnames(richness.rare) <- c('Chao1','se.chao1', 'Shannon', 'Overwintered')
  richness.rare$Labels <- rownames(richness.rare)
  ps.prop <- transform_sample_counts(ps, function(x) x/sum(x)*100)
  
  ad.test.df <- richness.rare[,c('Shannon', 'Chao1')]
  ad.test.df <- cbind(ad.test.df,
                      sample_data(ps.prop))
  
  ad.wilcox.shannon <- ad.test.df %>%
    wilcox_test(Shannon ~ Overwintered) %>%
    adjust_pvalue() %>%
    add_significance() %>%
    filter(p.adj.signif != "ns")
  
  ad.wilcon.Chao1 <- ad.test.df %>%
    wilcox_test(Chao1 ~ Overwintered) %>%
    adjust_pvalue() %>%
    add_significance() %>%
    filter(p.adj.signif != "ns")
  ad.test.df$Overwintered = factor(ad.test.df$Overwintered)
  
  print(ad.test.df %>%
          wilcox_test(Shannon ~ Overwintered) %>%
          adjust_pvalue() %>%
          add_significance())
  print(ad.test.df %>%
          wilcox_test(Chao1 ~ Overwintered) %>%
          adjust_pvalue() %>%
          add_significance())
  shannon.plot.sig <- ggboxplot(ad.test.df, x = "Overwintered", y = "Shannon",
                                color = "black", palette = primary_color_list, fill = "Overwintered",breaks = c('No', 'Yes')) +
    stat_pvalue_manual(ad.wilcox.shannon, label = "p.adj.signif", y.position = (max(ad.test.df$Shannon))+0.1, step.increase = 0.1) + 
    ylab('Shannon Index') +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.x = element_blank()) + 
    grids(linetype = "dashed")
  
  Chao1.plot.sig <- ggboxplot(ad.test.df, x = "Overwintered", y = "Chao1",
                              color = "black", palette = primary_color_list, fill = "Overwintered", breaks = c('No', 'Yes')) + 
    stat_pvalue_manual(ad.wilcon.Chao1, label = "p.adj.signif", y.position = (max(ad.test.df$Chao1)+1), step.increase = 0.1) + 
    ylab('Chao1 Index') +
    theme(legend.position = 'right',
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.x = element_blank()) + 
    grids(linetype = "dashed")
  
  diversity.plots.sig <- plot_grid(shannon.plot.sig, Chao1.plot.sig, 
                                   labels = c('A', 'B'), align = 'h',
                                   rel_widths = c(1.5, 2))  
  png(paste('./ov_',clade,'_',loc,'_alpha_div.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(diversity.plots.sig)
  dev.off()
  return(list(shannon.plot.sig,Chao1.plot.sig))
}
beta_diversity_ov <- function(ps,clade="None",loc="",prev_th=0.1,type_dist='bray'){
  if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
  # Prevalence filtering to 10% (default, specify if different) prevalence
  prevThreshold <- nsamples(ps) * prev_th
  
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
  
  # Generate nmds ordinations (Bray-Curtis dissimilarity)
  bc.nmds <- metaMDS(ord.asvs,disatnce = type_dist,
                     autotransform = F,
                     trymax = 1000)
  bc.nmds
  
  # Shepards test/goodness of fit
  gn <- goodness(bc.nmds) # Produces a results of test statistics for goodness of fit for each point
  png(paste('./ov_',clade,'_',loc,'_beta_div_stressplot.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(stressplot(bc.nmds))
  dev.off()
  
  # Calculate distance matrices
  bc.dist <- metaMDSredist(bc.nmds)
  bc.envfit.nmds <- envfit(bc.nmds ~ Overwintered,
                           data = ord.meta,
                           display = 'sites',
                           na.rm = T)
  
  # Extract fit data
  # NMDS data
  # Use ggvegan to extract the nmds data for plotting as ggplot object
  bc.envfit.nmds.fort <- fortify(bc.envfit.nmds)
  centroids.nmds <- subset(bc.envfit.nmds.fort, Type == 'Centroid')
  centroids.nmds$Label <- c('No', 'Yes')
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
                           Overwintered = ord.meta$Overwintered)
  
  plot.new()
  bc.nmds.ellipse <- ordiellipse(bc.nmds, 
                                 ord.meta$Overwintered, 
                                 kind = "sd", 
                                 conf = 0.95,
                                 object = T)
  bc.nmds.ellipse.df <- data.frame()
  ord.meta$Overwintered <- as.factor(ord.meta$Overwintered)
  for(tp in levels(ord.meta$Overwintered)){
    bc.nmds.ellipse.df <- rbind(bc.nmds.ellipse.df, 
                                cbind(as.data.frame(with(bc.nmds.df[bc.nmds.df$Overwintered==tp,],
                                                         veganCovEllipse(bc.nmds.ellipse[[tp]]$cov,
                                                                         bc.nmds.ellipse[[tp]]$center,
                                                                         bc.nmds.ellipse[[tp]]$scale)))
                                      ,Overwintered=tp))
  }
  
  
  
  ord.bray.nmds.plot <- plot_ordination(
    ps2.comp,
    bc.nmds,
    color = 'Overwintered', shape = 'Overwintered') +
    theme_cowplot() + 
    background_grid(major = 'xy', minor = 'none') +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = primary_color_list)
  
  bc.plot <- plot(ord.bray.nmds.plot +
    geom_path(data = bc.nmds.ellipse.df,
              aes(x = NMDS1, y = NMDS2, color = Overwintered),
              size = 1, 
              linetype = 2,
              inherit.aes = F))
  
  bc.adonis <- adonis2(bc.dist ~ Overwintered,
                       data = ord.meta)
  
  png(paste('./ov_',clade,'_',loc,'_beta_div.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(bc.adonis)
  dev.off()
  
  bc.disper <- betadisper(bc.dist, ord.meta$Overwintered)
  print(bc.adonis)
  pertest <- permutest(bc.disper,permutations = 999)
  print(pertest)
  bc.pw <- pairwise.adonis(
    otu_table(ps),
    sample_data(ps)$Overwintered,
    sim.method = "bray",
    p.adjust.m = "none",
    
  )
  print(bc.pw)
  capture.output(bc.disper, file = paste('./ov_',clade,'_',loc,'_beta_div_bcdisper.txt',sep=''))
  capture.output(pertest, file = paste('./ov_',clade,'_',loc,'_beta_div_permutest.txt',sep=''))
  capture.output(bc.pw, file = paste('./ov_',clade,'_',loc,'_beta_div_pwadonis.txt',sep=''))
  return(bc.plot)
}
beta_diversity_pop <- function(ps,clade="None",loc="",prev_th=0.1,type_dist='bray', seq=""){
  if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
  # Prevalence filtering to 10% (default, specify if different) prevalence
  prevThreshold <- nsamples(ps) * prev_th
  
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
  
  # Generate nmds ordinations (Bray-Curtis dissimilarity)
  bc.nmds <- metaMDS(ord.asvs,disatnce = type_dist,
                     autotransform = F,
                     trymax = 1000)
  bc.nmds
  
  # Shepards test/goodness of fit
  gn <- goodness(bc.nmds) # Produces a results of test statistics for goodness of fit for each point
  png(paste('./pop_',clade,'_',loc,'_beta_div_stressplot.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(stressplot(bc.nmds))
  dev.off()
  
  # Calculate distance matrices
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
                           Population = factor(ord.meta$Population, 
                                               levels=c('BS','TN','SM','DN')))
  
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
  
  
  
  ord.bray.nmds.plot <- plot_ordination(
    ps2.comp,
    bc.nmds,
    color = 'Population', shape = 'Population') +
    theme_cowplot() + 
    background_grid(major = 'xy', minor = 'none') +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = primary_color_list,breaks = c('BS','TN','SM','DN'))
  
  bc.plot <- ord.bray.nmds.plot +
    geom_path(data = bc.nmds.ellipse.df,
              aes(x = NMDS1, y = NMDS2, color = Population),
              size = 1, 
              linetype = 2,
              inherit.aes = F) + 
    grids(linetype = "dashed") 
  
  bc.adonis <- adonis2(bc.dist ~ Population,
                       data = ord.meta)
  
  png(paste('./pop_',clade,'_',loc,'_beta_div.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(bc.plot)
  dev.off()
  
  bc.disper <- betadisper(bc.dist, ord.meta$Population)
  print(bc.adonis)
  pertest <- permutest(bc.disper,permutations = 999)
  print(pertest)
  bc.pw <- pairwise.adonis(
    otu_table(ps),
    sample_data(ps)$Population,
    sim.method = "bray",
    p.adjust.m = "none",
    
  )
  print(bc.pw)
  capture.output(bc.disper, file = paste('./pop_',clade,'_',loc,'_beta_div_bcdisper.txt',sep=''))
  capture.output(pertest, file = paste('./pop_',clade,'_',loc,'_beta_div_permutest.txt',sep=''))
  capture.output(bc.pw, file = paste('./pop_',clade,'_',loc,'_beta_div_pwadonis.txt',sep=''))
  return(bc.plot)
}
DESeq_ov.pop <- function(ps,pop,print.all=F,clade="None"){
  # multi factorial analysis
  if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
  rownames(sample_data(ps)) <- sample_data(ps)$Subject
  dds <- phyloseq_to_deseq2(ps, ~ Population + Overwintered)
  keep <- rowSums(counts(dds)) >= 10
  ddsMF <- dds[keep,]
  levels(ddsMF$Overwintered)
  levels(ddsMF$Population)
  ddsMF$group <- factor(paste0(ddsMF$Overwintered, ddsMF$Population))
  levels(ddsMF$group)
  design(ddsMF) <- formula(~ group)
  ddsMF <- DESeq(ddsMF, test="Wald", fitType="parametric")
  resMF <- results(ddsMF)
  resultsNames(ddsMF)
  resMFType <- results(ddsMF, contrast=c("group", paste("Yes",pop,sep=''),
                                         paste("No",pop,sep='')))
  vsd <- varianceStabilizingTransformation(ddsMF, blind=FALSE) # vst normalization
  
  resMFType <- resMFType[order(resMFType$padj),] # order by adjusted p-value
  alpha = 0.05
  sigtab = resMFType[which(resMFType$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  dim(sigtab)
  sigtab <- sigtab[-14] # removing doubled species column

  tax.sig <- as.data.frame(cbind(rownames(sigtab),tax_table(ps)[rownames(sigtab),]))[,-8] # get tax info and removing subspecies
  
  if (print.all==T){
    for (i in 1:dim(tax.sig)[1]){ # plot all significant ASVs in relation to the samples as barplots
      if (is.na(tax.sig[i,]$Family)){titolo <- tax.sig[i,]$Order}
      else {
        if (is.na(tax.sig[i,]$Genus)){titolo <- tax.sig[i,]$Family}
        else {
          if (is.na(tax.sig[i,]$Species)){titolo <- tax.sig[i,]$Genus}
          else {titolo <- paste(tax.sig[i,]$Genus,tax.sig[i,]$Species)}
        }
      }
      png(paste(titolo,pop,'.png',sep='_'), width=6.5,units="in", height=5, res=1200)
      print(barplot(assay(ddsMF)[c(rownames(sigtab[i,])),c(rownames(df))],las=2,main=titolo))
      dev.off()
    } 
  }
  
  # log2fold change plot
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  write_xlsx(cbind(" "=rownames(sigtab),sigtab),paste(pop,'_ov_sigtaxa.xlsx',sep=''))
  png(paste('log2fold_',pop,'_DESeq2.png',sep=''), width=6.5,units="in", height=5, res=1200) # uncomment if you want to save the plot
  print(ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_jitter(size=3, width = 0.2) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)))
  dev.off() # uncomment if you want to save the plot
  
  # presence/absence plot of the most differential present ASVs 
  # prova con nome asv al posto delle asv
  # cambiare colore in verde e arancio
  df <- as.data.frame(colData(ddsMF)[,c("Overwintered","Population")])
  df$Overwintered <- factor(df$Overwintered, levels=c('No','Yes'))
  df <- df[order(factor(df$Overwintered, levels=c('No','Yes'))),]
  df <- df[which(df$Population == pop),]
  df <- df %>% dplyr::select(Overwintered)
  
  sig_rows <- sigtab %>% 
    mutate(genus = ifelse(!is.na(Species), paste(as.character(Genus),as.character(Species)), as.character(Genus))) %>% 
    mutate(genus = ifelse(!is.na(genus), genus, as.character(Family))) %>% 
    mutate(genus = ifelse(!is.na(genus), genus, as.character(Order))) %>%
    mutate(genus = ifelse(!is.na(genus), genus, as.character(Class))) %>%
    select(genus)
  sig_rows$genus <- ifelse(grepl("Rhizobium",sig_rows$genus),"Rhizobium",sig_rows$genus)
  rownames(sig_rows) <- rownames(assay(vsd)[c(rownames(sigtab)),c(rownames(df))])
  my_colors <- c("No" = "#1B9E779A", "Yes" = "#D95F039A")
  plot_deseq <- pheatmap(assay(vsd)[c(rownames(sigtab)),c(rownames(df))], cluster_rows=FALSE, show_rownames=T,
                         cluster_cols=FALSE, annotation_col=df, annotation_colors=list(Overwintered=my_colors),
                         labels_row=sig_rows$genus)
  png(paste(pop,'_DESeq2_sig.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(plot_deseq)
  dev.off()
  
  
  # uncomment the following three lines in case you are interested in looking at 
  # the ASV counts instead of the significant differential present ASVs
  
  # select <- order(rowMeans(counts(ddsMF,normalized=TRUE)), decreasing=TRUE)[1:20]
  # pheatmap(assay(vsd)[c(rownames(resMFType[select,])),], cluster_rows=FALSE, show_rownames=FALSE,
  #         cluster_cols=FALSE, annotation_col=df)
  
  return(list(sigtab,plot_deseq))
}
picrust2_analysis <- function(ps,p2EC,p2KO,p2PW,var,pop='full_data'){
  # Subset from the ps
  p2EC = p2EC[,sample_names(ps)]
  p2KO = p2KO[,sample_names(ps)]
  p2PW = p2PW[,sample_names(ps)]
  
  # EC
  set.seed(12345)
  system.time({
    aldex2_EC = aldex(p2EC, sample_data(ps)[[var]], mc.samples = 500, test = "t", 
                          effect = TRUE, denom = "iqlr", verbose = TRUE)
  })
  # KO
  set.seed(12345)
  system.time({
    aldex2_KO = aldex(p2KO, sample_data(ps)[[var]], mc.samples = 500, test = "t", 
                          effect = TRUE, denom = "iqlr", verbose = TRUE)
  })
  # Pathway
  # KO
  set.seed(12345)
  system.time({
    aldex2_PW = aldex(p2PW, sample_data(ps)[[var]], mc.samples = 500, test = "t", 
                          effect = TRUE, denom = "iqlr", verbose = TRUE)
  })
  
  png(paste("./ALDEx2_picrust2_effect_",pop,".png",sep=''), width = 6, height = 6, units = "in", res = 300)
  par(mfrow = c(2,2))
  hist(aldex2_EC$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
  hist(aldex2_KO$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
  hist(aldex2_PW$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
  invisible(dev.off())
  
  png(paste("./ALDEx2_picrust2_MW_MA_",pop,".png",sep=''), width = 6, height = 8, units = "in", res = 300)
  par(mfrow = c(3,2))
  aldex.plot(aldex2_EC, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
             called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
  title(main = "(EC) MW Plot")
  
  aldex.plot(aldex2_EC, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
             called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference")
  title(main = "(EC) MA Plot")
  
  aldex.plot(aldex2_KO, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
             called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
  title(main = "(KO) MW Plot")
  
  aldex.plot(aldex2_KO, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
             called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
  title(main = "(KO) MA Plot")
  
  aldex.plot(aldex2_PW, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
             called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
  title(main = "(PW) MW Plot")
  
  aldex.plot(aldex2_PW, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
             called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
  title(main = "(PW) MA Plot")
  invisible(dev.off())
  
  png(paste("./ALDEx2_picrust2_P_adjP_",pop,".png",sep=''), width = 6, height = 8, units = "in", res = 300)
  par(mfrow = c(3,2))
  plot(aldex2_EC$effect, aldex2_EC$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
       xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
  points(aldex2_EC$effect, aldex2_EC$we.eBH, cex = 0.5, col = "red", pch = 19)
  abline(h = 0.05, lty = 2, col = "grey")
  legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))
  
  plot(aldex2_EC$diff.btw, aldex2_EC$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
       xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
  points(aldex2_EC$diff.btw, aldex2_EC$we.eBH, cex = 0.5, col = "red", pch = 19)
  abline(h = 0.05, lty = 2, col = "grey")
  
  plot(aldex2_KO$effect, aldex2_KO$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
       xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
  points(aldex2_KO$effect, aldex2_KO$we.eBH, cex = 0.5, col = "red", pch = 19)
  abline(h = 0.05, lty = 2, col = "grey")
  legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))
  
  plot(aldex2_KO$diff.btw, aldex2_KO$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
       xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
  points(aldex2_KO$diff.btw, aldex2_KO$we.eBH, cex = 0.5, col = "red", pch = 19)
  abline(h = 0.05, lty = 2, col = "grey")
  
  plot(aldex2_PW$effect, aldex2_PW$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
       xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
  points(aldex2_PW$effect, aldex2_PW$we.eBH, cex = 0.5, col = "red", pch = 19)
  abline(h = 0.05, lty = 2, col = "grey")
  legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))
  
  plot(aldex2_PW$diff.btw, aldex2_PW$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
       xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
  points(aldex2_PW$diff.btw, aldex2_PW$we.eBH, cex = 0.5, col = "red", pch = 19)
  abline(h = 0.05, lty = 2, col = "grey")
  invisible(dev.off())
  
  df_EC <- aldex2_EC %>% tibble::rownames_to_column(var = "EC") %>% 
    inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)
  
  df_KO = aldex2_KO %>% tibble::rownames_to_column(var = "KO") %>% 
    inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)
  
  df_PW = aldex2_PW %>% tibble::rownames_to_column(var = "Pathway") %>% 
    inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)
  
  write.table(df_EC, file = paste("./ALDEx2_picrust2_EC_results_",pop,".tsv",sep=''), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  write.table(df_KO, file = paste("./ALDEx2_picrust2_KO_results_",pop,".tsv",sep=''), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  write.table(df_PW, file = paste("./ALDEx2_picrust2_Pathway_results_",pop,".tsv",sep=''), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  return()
}
pantoea_check <- function(ps){
  ps2.bar <- ps
  # Merged phyloseq plot (100% total abundance)
  ps2.bar.merge <- merge_samples(ps2.bar, group = 'Overwintered')
  ## Warning in asMethod(object): NAs introduced by coercion
  sample_data(ps2.bar.merge)$Overwintered <- sample_names(ps2.bar.merge)
  ps2.bar.merge.comp <- transform_sample_counts(ps2.bar.merge,  
                                                function(x) x/sum(x) * 100)
  ps2.bp.Genus <- tax_glom(ps2.bar.merge.comp, 'Genus', NArm = F)
  pantoea <- subset_taxa(ps2.bp.Genus, Genus=="Pantoea")
  df <- cbind(sample_data(pantoea),otu_table(pantoea))
  colnames(df)[7] <- '% Pantoea'
  df$`% Pantoea` <- round(df$`% Pantoea`, digits=2)
  bp<-ggplot(data=df, aes(x=Overwintered, y=`% Pantoea`,fill=Overwintered)) +
    geom_bar(stat="identity", color="black")+
    geom_text(aes(label=`% Pantoea`), vjust=-0.3, size=3.5)+
    theme_minimal()
  png('./ov_pantoea_perc_barplot.png', width=6.5,units="in", height=5, res=1200)
  print(bp)
  dev.off()
  ps.bp.Genus <- tax_glom(ps, 'Genus', NArm = F)
  wol_puglia <- subset_taxa(ps.bp.Genus, Genus=="Pantoea")
  df <- cbind(sample_data(wol_puglia),otu_table(wol_puglia))
  colnames(df)[7] <- 'Pantoea'
  
  bp <- ggplot(data=df, aes(x=Subject, y=Pantoea,fill=Overwintered)) +
    geom_bar(stat="identity")+
    geom_text(aes(label=Pantoea), vjust=-0.3, size=3.5)+
    theme_minimal()
  png('./ov_pantoea_barplot.png', width=6.5,units="in", height=5, res=1200)
  print(bp)
  dev.off()
  print(wilcox.test(df$Pantoea~df$Overwintered))
}
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni'){
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[grepl(paste(c(co[1, elem], co[2, elem]), collapse = "|"), factors),],method=sim.method)}
    
    ad = adonis2(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$F[1]);
    R2 = c(R2,ad$R2[1]);
    p.value = c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 
beta_diversity_pop_general <- function(ps,clade="None",loc="",prev_th=0.1,type_dist='bray'){
  if (clade != "None"){ps <- tax_glom(ps,clade,NArm=F)}
  # Prevalence filtering to 10% (default, specify if different) prevalence
  prevThreshold <- nsamples(ps) * prev_th
  l <- levels(sample_data(ps)$Population)
  
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
  
  # Generate nmds ordinations (Bray-Curtis dissimilarity)
  bc.nmds <- metaMDS(ord.asvs,disatnce = type_dist,
                     autotransform = F,
                     trymax = 1000)
  bc.nmds
  
  # Shepards test/goodness of fit
  gn <- goodness(bc.nmds) # Produces a results of test statistics for goodness of fit for each point
  png(paste('./pop_',clade,'_',loc,'_beta_div_stressplot.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(stressplot(bc.nmds))
  dev.off()
  
  # Calculate distance matrices
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
  centroids.nmds$Label <- l
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
                           Population = factor(ord.meta$Population, 
                                               levels=l))
  
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
  
  
  
  ord.bray.nmds.plot <- plot_ordination(
    ps2.comp,
    bc.nmds,
    color = 'Population', shape = 'Population') +
    theme_cowplot() + 
    background_grid(major = 'xy', minor = 'none') +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = primary_color_list,breaks = l)
  
  bc.plot <- ord.bray.nmds.plot +
    geom_path(data = bc.nmds.ellipse.df,
              aes(x = NMDS1, y = NMDS2, color = Population),
              size = 1, 
              linetype = 2,
              inherit.aes = F) + 
    grids(linetype = "dashed") 
  
  bc.adonis <- adonis2(bc.dist ~ Population,
                       data = ord.meta)
  
  png(paste('./pop_',clade,'_',loc,'_beta_div.png',sep=''), width=6.5,units="in", height=5, res=1200)
  print(bc.plot)
  dev.off()
  
  bc.disper <- betadisper(bc.dist, ord.meta$Population)
  print(bc.adonis)
  pertest <- permutest(bc.disper,permutations = 999)
  print(pertest)
  bc.pw <- pairwise.adonis(
    otu_table(ps),
    sample_data(ps)$Population,
    sim.method = "bray",
    p.adjust.m = "none",
    
  )
  print(bc.pw)
  capture.output(bc.disper, file = paste('./pop_',clade,'_',loc,'_beta_div_bcdisper.txt',sep=''))
  capture.output(pertest, file = paste('./pop_',clade,'_',loc,'_beta_div_permutest.txt',sep=''))
  capture.output(bc.pw, file = paste('./pop_',clade,'_',loc,'_beta_div_pwadonis.txt',sep=''))
  return(bc.plot)
}