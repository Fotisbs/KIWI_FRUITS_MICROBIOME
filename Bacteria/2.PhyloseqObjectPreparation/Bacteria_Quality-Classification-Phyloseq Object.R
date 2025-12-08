##Quality Control, classification and phyloseq object construction##

# install the necessary packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")


#### package loading ----
# load the package and check the package version
library(dada2); packageVersion("dada2")

mymultthread <- 56


##input file path and storing in variables##
##change to the directory containing the fastq files after unzipping##
path1 <- ("~/ex.Bacteria")

# lists files in the path
list.files(c(path1))

# set the variables containing all the forward and the reverse paths to the files of interest with the list.files command
fnFs <- sort(list.files(c(path1), pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(c(path1), pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub("_lib_.+","",basename(fnFs))

# Plot the per-base qualities of all the samples in a single file
pdf("Initial Quality Files.pdf",onefile = T)
for (i in c(1:length(fnFs))) {
  print(i)
  plot1 <- plotQualityProfile(c(fnFs[i], fnRs[i]))
  print(plot1)
}
dev.off()


##Sequence quality filtering and control, error modelling, and dereplication##
##Set the file paths where the quality controlled sequences are going to be saved##
filtFs <- file.path("filtered2", paste(sample.names, "_F_filt.fastq.gz", sep = ""))
filtRs <- file.path("filtered2", paste(sample.names, "_R_filt.fastq.gz", sep = ""))

##filter the sequences and save then in the folders provided above and get their statistics in a table##
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft = 11, trimRight = 50, maxN=0,  
                     maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, 
                     multithread=TRUE, matchIDs=TRUE)
head(out)

##learns error rates using a machine learning algorithm##
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize the estimated error rates in a single pdf
pdf("Estimated Error Rates.pdf",onefile = T)
plotErrF <- plotErrors(errF, nominalQ=TRUE)
print(plotErrF)
plotErrR <- plotErrors(errR, nominalQ=TRUE)
print(plotErrR)
dev.off()

##Dereplication of each one of the red pairs to unique sequences (collapsing of the identical sequences for each pair per sample)##
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

##Name the derep-class objects by the sample names##
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Sample composition inference after read correction##
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

##Merge read pairs retaining the per sequence sample information##
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

##Construct the sequence table, remove the chimeras, and create a summary##
##Construct sequence table (ASV table)##
seqtab <- makeSequenceTable(mergers)

View(seqtab)

#Number of samples
nrow(seqtab)
#Number of sequence Variants-ASVs##
ncol(seqtab)

#Distribution of sequence lengths##
table(nchar(getSequences(seqtab)))

##Chimera removal (consensus)##
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mymultthread, verbose=TRUE)

##View(seqtab.nochim)##
dim(seqtab.nochim)

##Record the portion of good sequences out of the total prior the chimera removal##
sum(seqtab.nochim)/sum(seqtab)

##Track reads through the pipeline##
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

head(track)

##Save the read quality control data##
write.table(track, file = "readQC_Bacteria.txt", col.names = NA, sep = "\t", quote = FALSE)

##Taxonomically classify the sequences##
taxa <- assignTaxonomy(seqtab.nochim, minBoot = 80, "/mnt/4abcb8af-e2eb-47e2-86d6-eed71f6d5304/00_home_guests_do_not_move/fotisbs/wshp_foldrs/dbs/old_dep/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)

## add species assignment##
taxa2 <- addSpecies(taxa, "/mnt/4abcb8af-e2eb-47e2-86d6-eed71f6d5304/00_home_guests_do_not_move/fotisbs/wshp_foldrs/dbs/old_dep/silva_species_assignment_v138.fa.gz")

write.table(taxa, file = "taxa.txt", col.names = NA, sep = "\t", quote = F)
write.table(taxa2, file = "taxa2.txt", col.names = NA, sep = "\t", quote = F)

##Proceed with analysis, load the phyloseq package##
library(phyloseq); packageVersion("phyloseq")

library(ggplot2); packageVersion("ggplot2")

#The following globally sets the theme of the plots##
theme_set(theme_bw())

##Load the experimental set up file (samdf)##
samdfBacteria <- read.table(file = "samdf.txt", header = TRUE, row.names = 1, sep = "\t")

View(samdfBacteria)

##Construct the phyloseq object## 
Bacteria_KIWI_2024 <- phyloseq(otu_table(seqtab.nochim, 
                                            taxa_are_rows=FALSE),
                                  sample_data(samdfBacteria),
                                  tax_table(taxa2))

##Replace the taxon names (the sequences of each ASV with something easier to read and save the sequence info in the object)##
##load the Biostrings and stringr packages##
library("Biostrings")

sequences <- Biostrings::DNAStringSet(taxa_names(Bacteria_KIWI_2024))

names(sequences) <- taxa_names(Bacteria_KIWI_2024)

ps <- merge_phyloseq(Bacteria_KIWI_2024, sequences)

Bacteria_KIWI_2024 <- ps

##The following command (uses the stringr::str_pad for padding leading zeros) performs the actual renaming##
library(stringr)

taxa_names(Bacteria_KIWI_2024) <- paste("ASV",str_pad(1:length(taxa_names(Bacteria_KIWI_2024)),4, pad = "0"),sep = "")

rank_names(Bacteria_KIWI_2024)

Bacteria_KIWI_2024

##SAVE THE FINAL PHYLOSEQ OBJECT##
saveRDS(Bacteria_KIWI_2024, file = "Bacteria_KIWI_2024.RDS")

##Check the object features##
View(data.frame(otu_table(Bacteria_KIWI_2024)))
View(data.frame(tax_table(Bacteria_KIWI_2024)))
View(data.frame(sample_data(Bacteria_KIWI_2024)))

##Create tables,and check number of features for each taxa##
table(tax_table(Bacteria_KIWI_2024)[, "Kingdom"], exclude = NULL)
table(tax_table(Bacteria_KIWI_2024)[, "Phylum"], exclude = NULL) 
table(tax_table(Bacteria_KIWI_2024)[, "Class"], exclude = NULL) 
table(tax_table(Bacteria_KIWI_2024)[, "Order"], exclude = NULL)
table(tax_table(Bacteria_KIWI_2024)[, "Family"], exclude = NULL)
table(tax_table(Bacteria_KIWI_2024)[, "Genus"], exclude = NULL)
table(tax_table(Bacteria_KIWI_2024)[, "Species"], exclude = NULL)


##Further, features with ambiguous annotation, low abundance or non target taxa removed from phyloseq object##
##For our analysis removed Eukaryota, Archaea, NA##
Bacteria_KIWI_2024_Clean <- subset_taxa(Bacteria_KIWI_2024, !is.na(Kingdom) & !Kingdom %in% (c("", "uncharacterized", "Eukaryota","Archaea")))

Bacteria_KIWI_2024_Clean <- prune_taxa(taxa_sums(Bacteria_KIWI_2024_Clean)>0, Bacteria_KIWI_2024_Clean)

##Also remove mitochondria and chloroplasts##
Bacteria_KIWI_2024_Clean <- subset_taxa(Bacteria_KIWI_2024_Clean, !Order %in% c("Chloroplast"))

Bacteria_KIWI_2024_Clean <- subset_taxa(Bacteria_KIWI_2024_Clean, !Family %in% c("Mitochondria"))

Bacteria_KIWI_2024_Clean <- prune_taxa(taxa_sums(Bacteria_KIWI_2024_Clean)>0, Bacteria_KIWI_2024_Clean)

##Replace the NA names with the classified leftmost (higher level annotated) taxa##
Bacteria_KIWI_2024_Annotated <- Bacteria_KIWI_2024_Clean

for(i in 1:nrow(tax_table(Bacteria_KIWI_2024_Annotated))){
  for(j in 2:ncol(tax_table(Bacteria_KIWI_2024_Annotated))){
    if(is.na(tax_table(Bacteria_KIWI_2024_Annotated)[i,j])){
      tax_table(Bacteria_KIWI_2024_Annotated)[i,j] <- tax_table(Bacteria_KIWI_2024_Annotated)[i,j-1]
    }
  }
}

####SAVE THE FINAL Annotated phyloseq object and continue with the statistical analysis scripts##
saveRDS(Bacteria_KIWI_2024_Annotated, file = "Bacteria_KIWI_2024_Annotated.RDS")