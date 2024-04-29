##### DADA2 Tutorial (1.16) https://benjjneb.github.io/dada2/tutorial.html #####

# rm(list=ls(all=TRUE))
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(DESeq2); packageVersion("DESeq2")
library(knitr)

# Set the path where to find Illumina MiSeq fastq files 
setwd('D:/PhD/PhD.AES.Dsuz.Riccardo.Piccinno/H_halys/metaBUG/Sequenze_MiSeq')
path <- "./FEM_AMPLICON_RUN178"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) # Get all forward reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) # Get all reverse reads
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Considerations for your own data: The string manipulations may have to be modified if your filename format is different

# Inspect read quality profiles
plotQualityProfile(fnFs[1:2]) # Quality of 1st and 2nd forward reads 
# The green line shows the mean quality score at each position
# The orange line shows the quartiles of the quality score distribution
# The red lines shows the proportion of reads that extend at least that position 
# The grey squares show a heat map of the frequency of each quality score at each base position
plotQualityProfile(fnRs[1:2]) # Quality of 1st and 2nd reverse reads 
# We can see the quality of reverse reads are significantly worse quality at the end
# However, DADA2 is robust to that kind of situation, since it implements this 
# information in the error model. 
# In this case, reverse reads will be truncated at 160th position.
# When using V3-V4 amplicons, 'truncLen' (see next step) must be larger than large 
# enough to maintain '20 + biological.length.variation nucleotides' of overlap between them.

# Filter and Trim
filtFs <- file.path('.', "FEM_AMPLICON_RUN178_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path('.', "FEM_AMPLICON_RUN178_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# The following step takes approximately 1 minute
# Sys.time()
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
# On Windows set multithread=FALSE, on linux set multithread=TRUE
# Sys.time()
head(out)
# maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE 
# parameter sets the maximum number of “expected errors” allowed in a read, 
# which is a better filter than simply averaging quality scores.
# If you want to speed up downstream computation, consider tightening maxEE. 
# If too few reads are passing the filter, consider relaxing maxEE, perhaps 
# especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen 
# to remove low quality tails. Remember though, when choosing truncLen for 
# paired-end reads you must maintain overlap after truncation in order to merge them later.
# In case of ITS follow the proper guide.


# Learn error rates
# In this step we will build an error model from the data.
# The following steps take approximately 3 minutes (33514080 bases in 139642 reads 
# from 20 samples repeated for the forward and the reverse)
# Sys.time()
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# Sys.time()
plotErrors(errF, nominalQ=TRUE)
# Plot showing error rates of all possible nucleotide trnasitions.
# The black line shows the estimated error rates after convergence of the 
# machine-learning algorithm. The red line shows the error rates expected under 
# the nominal definition of the Q-score.
# Everything looks reasonable in this case.

# Sample inference
# The following steps take approximately 1 minute
# Sys.time()
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# Sys.time()
# Need to do some minor changes in case of IonTorrent or 454 pyrosequencing
dadaFs[[1]]
# 128 true sequence variants from the 1979 unique sequences in the first sample 
# inferred by the algorithm.
# help("dada-class") for some info about the output

# Merge paired ends
# Merging is performed by aligning the denoised forward reads with the reverse-complement 
# of the corresponding denoised reverse reads, and then constructing the merged “contig” 
# sequences. By default, merged sequences are only output if the forward and reverse reads 
# overlap by at least 12 bases, and are identical to each other in the overlap region 
# (but these conditions can be changed via function arguments).
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# The mergers object is a list of data.frames from each sample. Each data.frame 
# contains the merged $sequence, its $abundance, and the indices of the $forward 
# and $reverse sequence variants that were merged. Paired reads that did not 
# exactly overlap were removed by mergePairs, further reducing spurious output.
# Most of your reads should successfully merge. If that is not the case upstream 
# parameters may need to be revisited: Did you trim away the overlap between your reads?
# Non-overlapping reads are supported, but not recommended, with mergePairs(..., justConcatenate=TRUE).

# Construct seqeunce table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# The sequence table is a matrix with rows corresponding to (and named by) the samples, 
# and columns corresponding to (and named by) the sequence variants. This table 
# contains 293 ASVs, and the lengths of our merged sequences all fall within the 
# expected range for this V4 amplicon.
# Sequences that are much longer or shorter than expected may be the result of 
# non-specific priming. You can remove non-target-length sequences from your sequence 
# table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). 
# This is analogous to “cutting a band” in-silico to get amplicons of the targeted length.

# Remove chimeras
# Chimeric sequences are identified if they can be exactly reconstructed by combining 
# a left-segment and a right-segment from two more abundant “parent” sequences.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# 4% of chimeras abundance
# Most of your reads should remain after chimera removal (it is not uncommon for
# a majority of sequence variants to be removed though). If most of your reads 
# were removed as chimeric, upstream processing may need to be revisited. In almost 
# all cases this is caused by primer sequences with ambiguous nucleotides that 
# were not removed prior to beginning the DADA2 pipeline.

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
# The majority of reads were kept!
# There should no step in which a majority of reads are lost. 
# If a majority of reads failed to merge, you may need to revisit the 'truncLen' 
# parameter used in the filtering step and make sure that the truncated reads span 
# your amplicon. 
# If a majority of reads were removed as chimeric, you may need to revisit the 
# removal of primers, as the ambiguous nucleotides in unremoved primers interfere 
# with chimera identification.

# Assigning taxonomy
# Naive Bayes classifier method used.
# The assignTaxonomy function takes as input a set of sequences to be classified 
# and a training set of reference sequences with known taxonomy, and outputs taxonomic 
# assignments with at least minBoot bootstrap confidence.
# Formatted training fastas for the RDP training set, GreenGenes clustered at 97% 
# identity, and the Silva reference database are maintined.
# Before proceeding, download the silva_nr_v132_train_set.fa.gz file, and place 
# it in the directory with the fastq files.
taxa <- assignTaxonomy(seqtab.nochim, "./silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "./silva_species_assignment_v138.1.fa.gz")
# Exact matching (or 100% identity) is the only appropriate way to assign species 
# to 16S gene fragments. Currently, species-assignment training fastas are available 
# for the Silva and RDP 16S databases.
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# alternative with DECIPHER
# library(DECIPHER); packageVersion("DECIPHER")
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# Evaluate accuracy
#  To evaluate we use sequences of a mock community of 20 bacteria strains 
# unqs.mock <- seqtab.nochim["Mock",]
# unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
# mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
# match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
# cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
# cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

