##NMDS Triplot ggplot2 (Only Significant Metabolites included according to the highest R² values from envfit() (i.e., the best-fitting vectors)##

##Load the current libraries##
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

##In order to proceed the phyloseq object have to be merged according to replicates of the same treatments##  
Bacteria_KIWI_2024_Annotated

View(data.frame(sample_data(Bacteria_KIWI_2024_Annotated)))

Bacteria_KIWI_2024_An_Merged <- merge_samples(Bacteria_KIWI_2024_Annotated, "Tissue_MCP_Storage")

##Bacteria object also were tax glommed to Genus level for better results## 
Bacteria_KIWI_2024_An_Merged_Glommed <- tax_glom(Bacteria_KIWI_2024_An_Merged, taxrank = "Genus")

##After merging the sample data info are lost and you have to import again the characteristics of each merged treatment##
write.table(data.frame(sample_data(Bacteria_KIWI_2024_An_Merged_Glommed)), file="sample_data_Bacteria_KIWI_2024_An_Merged_Glommed.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("sample_data_Bacteria_KIWI_2024_An_Merged_Glommed.txt", header=T,sep = "\t",row.names = 1)

sample_data(Bacteria_KIWI_2024_An_Merged_Glommed) <- SampleDataNew67

View(data.frame(sample_data(Bacteria_KIWI_2024_An_Merged_Glommed)))

Bacteria_KIWI_2024_An_Merged_Glommed

##Transform to relative abundance %##
Bacteria_KIWI_2024_An_Merged_Glommed_100 <- transform_sample_counts(Bacteria_KIWI_2024_An_Merged_Glommed, function(OTU) 100*OTU/sum(OTU))

Bacteria_KIWI_2024_An_Merged_Glommed_100

##Ready to go with the NMDS Triplot Analysis, after uploading and scalling Metabolites tables (Amino_Acids, Organic_Acids, Carbohydrates, Phenolics)##

##Bacteria_KIWI vs Amino_Acids##
##Load and scale metabolites##
Amino_Acids_Triplot <- read.table("Amino_Acids_Triplot.txt", header = TRUE, sep = "\t", row.names = 1)
Amino_Acids_Triplot_scaled <- scale(Amino_Acids_Triplot)

##Phyloseq object##
physeq_raw <- Bacteria_KIWI_2024_An_Merged_Glommed_100

##Build OTU matrix##
otu_mat <- as(otu_table(physeq_raw), "matrix")
if (taxa_are_rows(physeq_raw)) {
  otu_mat <- t(otu_mat)
}

##NMDS ordination##
ord.nmds.bray <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)

##Metadata##
meta <- as(sample_data(physeq_raw), "data.frame")
meta$SampleID <- rownames(meta)

##Extract NMDS site scores##
site_scores <- as.data.frame(scores(ord.nmds.bray, display = "sites"))
site_scores$SampleID <- rownames(site_scores)
site_scores <- left_join(site_scores, meta, by = "SampleID")

##Fit metabolites to ordination##
fit_env <- envfit(ord.nmds.bray, Amino_Acids_Triplot_scaled)

##Extract top metabolites based on R²##
r2_values <- fit_env$vectors$r
top5_idx <- order(r2_values, decreasing = TRUE)[1:6]
top5_arrows <- fit_env$vectors$arrows[top5_idx, , drop = FALSE]
metabolite_arrows <- as.data.frame(top5_arrows * 0.3)
metabolite_arrows$Metabolite <- rownames(top5_arrows)

##Extract taxa scores##
species_scores <- scores(ord.nmds.bray, display = "species")
Top_25 <- names(sort(taxa_sums(physeq_raw), decreasing = TRUE)[1:20])
Top_25_in_otu <- intersect(Top_25, rownames(species_scores))
taxa_arrows <- as.data.frame(species_scores[Top_25_in_otu, ] * 1.8)
taxa_arrows$TaxonID <- rownames(taxa_arrows)

##Taxonomy labels##
tax_tbl <- as.data.frame(tax_table(physeq_raw))
tax_tbl_clean <- apply(tax_tbl, 2, function(x) gsub("^[a-zA-Z]__", "", x))
tax_tbl_clean <- as.data.frame(tax_tbl_clean)

taxa_arrows$Label <- sapply(taxa_arrows$TaxonID, function(tid) {
  family <- tax_tbl_clean[tid, "Family"]
  genus <- tax_tbl_clean[tid, "Genus"]
  order <- tax_tbl_clean[tid, "Order"]
  
  if (!is.na(family) && family != "" && !is.na(genus) && genus != "") {
    paste0(family, " ", genus)
  } else if (!is.na(family) && family != "") {
    paste0(family, " sp.")
  } else if (!is.na(order) && order != "") {
    paste0("Unclassified ", order)
  } else {
    "Unclassified"
  }
})

##ggplot2 colors and shapes##
myclos <- brewer.pal(n = 4, name = "RdBu")
shapes <- c(21, 22, 24, 25)

##Plot##
p <- ggplot(data = site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Cold_Storage, shape = Tissue), size = 4, color = "black") +
  scale_fill_manual(values = myclos) +
  scale_shape_manual(values = shapes[1:length(unique(site_scores$Tissue))]) +
  
## Top metabolite arrows##
  geom_segment(data = metabolite_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkblue", size = 1.2) +
  geom_text_repel(
    data = metabolite_arrows,
    aes(x = NMDS1, y = NMDS2, label = Metabolite),
    color = "darkblue", fontface = "bold", size = 4,
    max.overlaps = Inf) +
  
##Top taxa arrows##
  geom_segment(data = taxa_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "gray30", size = 0.8) +
  geom_text_repel(
    data = taxa_arrows,
    aes(x = NMDS1, y = NMDS2, label = Label),
    size = 3, fontface = "italic", color = "black",
    max.overlaps = Inf) +
  
  theme_minimal() +
  labs(title = "Bacteria–Amino Acids NMDS Triplot",
       subtitle = paste("NMDS Stress =", round(ord.nmds.bray$stress, 3)),
       fill = "Cold Storage",
       shape = "Tissue Type") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

##Save plot##
ggsave("Bacteria_Amino_Acids_NMDS_Triplot.pdf", plot = p, width = 8, height = 7)

##Bacteria_KIWI vs Organic_Acids##
# Load and scale metabolites##
Organic_Acids_Triplot <- read.table("Organic_Acids_Triplot.txt", header = TRUE, sep = "\t", row.names = 1)
Organic_Acids_Triplot_scaled <- scale(Organic_Acids_Triplot)

##Phyloseq object##
physeq_raw <- Bacteria_KIWI_2024_An_Merged_Glommed_100

##Build OTU matrix##
otu_mat <- as(otu_table(physeq_raw), "matrix")
if (taxa_are_rows(physeq_raw)) {
  otu_mat <- t(otu_mat)
}

##NMDS ordination##
ord.nmds.bray <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)

##Metadata##
meta <- as(sample_data(physeq_raw), "data.frame")
meta$SampleID <- rownames(meta)

##Extract NMDS site scores##
site_scores <- as.data.frame(scores(ord.nmds.bray, display = "sites"))
site_scores$SampleID <- rownames(site_scores)
site_scores <- left_join(site_scores, meta, by = "SampleID")

##Fit metabolites to ordination##
fit_env <- envfit(ord.nmds.bray, Organic_Acids_Triplot_scaled)

##Extract top metabolites based on R²##
r2_values <- fit_env$vectors$r
top5_idx <- order(r2_values, decreasing = TRUE)[1:9]
top5_arrows <- fit_env$vectors$arrows[top5_idx, , drop = FALSE]
metabolite_arrows <- as.data.frame(top5_arrows * 0.3)
metabolite_arrows$Metabolite <- rownames(top5_arrows)

##Extract taxa scores##
species_scores <- scores(ord.nmds.bray, display = "species")
Top_25 <- names(sort(taxa_sums(physeq_raw), decreasing = TRUE)[1:20])
Top_25_in_otu <- intersect(Top_25, rownames(species_scores))
taxa_arrows <- as.data.frame(species_scores[Top_25_in_otu, ] * 1.8)
taxa_arrows$TaxonID <- rownames(taxa_arrows)

##Taxonomy labels##
tax_tbl <- as.data.frame(tax_table(physeq_raw))
tax_tbl_clean <- apply(tax_tbl, 2, function(x) gsub("^[a-zA-Z]__", "", x))
tax_tbl_clean <- as.data.frame(tax_tbl_clean)

taxa_arrows$Label <- sapply(taxa_arrows$TaxonID, function(tid) {
  family <- tax_tbl_clean[tid, "Family"]
  genus <- tax_tbl_clean[tid, "Genus"]
  order <- tax_tbl_clean[tid, "Order"]
  
  if (!is.na(family) && family != "" && !is.na(genus) && genus != "") {
    paste0(family, " ", genus)
  } else if (!is.na(family) && family != "") {
    paste0(family, " sp.")
  } else if (!is.na(order) && order != "") {
    paste0("Unclassified ", order)
  } else {
    "Unclassified"
  }
})

##ggplot2 colors and shapes##
myclos <- brewer.pal(n = 4, name = "RdBu")
shapes <- c(21, 22, 24, 25)

##Plot##
p <- ggplot(data = site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Cold_Storage, shape = Tissue), size = 4, color = "black") +
  scale_fill_manual(values = myclos) +
  scale_shape_manual(values = shapes[1:length(unique(site_scores$Tissue))]) +
  
##Top metabolite arrows##
geom_segment(data = metabolite_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkblue", size = 1.2) +
  geom_text_repel(
    data = metabolite_arrows,
    aes(x = NMDS1, y = NMDS2, label = Metabolite),
    color = "darkblue", fontface = "bold", size = 4,
    max.overlaps = Inf) +
  
##Top taxa arrows##
geom_segment(data = taxa_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "gray30", size = 0.8) +
  geom_text_repel(
    data = taxa_arrows,
    aes(x = NMDS1, y = NMDS2, label = Label),
    size = 3, fontface = "italic", color = "black",
    max.overlaps = Inf) +
  
  theme_minimal() +
  labs(title = "Bacteria–Organic Acids NMDS Triplot",
       subtitle = paste("NMDS Stress =", round(ord.nmds.bray$stress, 3)),
       fill = "Cold Storage",
       shape = "Tissue Type") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

##Save plot##
ggsave("Bacteria_Organic_Acids_NMDS_Triplot.pdf", plot = p, width = 8, height = 7)


##Bacteria_KIWI vs Carbohydrates##
##Load and scale metabolites##
Carbohydrates_Triplot <- read.table("Carbohydrates_Triplot.txt", header = TRUE, sep = "\t", row.names = 1)
Carbohydrates_Triplot_scaled <- scale(Carbohydrates_Triplot)

##Phyloseq object##
physeq_raw <- Bacteria_KIWI_2024_An_Merged_Glommed_100

##Build OTU matrix##
otu_mat <- as(otu_table(physeq_raw), "matrix")
if (taxa_are_rows(physeq_raw)) {
  otu_mat <- t(otu_mat)
}

##NMDS ordination##
ord.nmds.bray <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)

##Metadata##
meta <- as(sample_data(physeq_raw), "data.frame")
meta$SampleID <- rownames(meta)

##Extract NMDS site scores##
site_scores <- as.data.frame(scores(ord.nmds.bray, display = "sites"))
site_scores$SampleID <- rownames(site_scores)
site_scores <- left_join(site_scores, meta, by = "SampleID")

##Fit metabolites to ordination##
fit_env <- envfit(ord.nmds.bray, Carbohydrates_Triplot_scaled)

##Extract top metabolites based on R²##
r2_values <- fit_env$vectors$r
top5_idx <- order(r2_values, decreasing = TRUE)[1:9]
top5_arrows <- fit_env$vectors$arrows[top5_idx, , drop = FALSE]
metabolite_arrows <- as.data.frame(top5_arrows * 0.3)
metabolite_arrows$Metabolite <- rownames(top5_arrows)

##Extract taxa scores##
species_scores <- scores(ord.nmds.bray, display = "species")
Top_25 <- names(sort(taxa_sums(physeq_raw), decreasing = TRUE)[1:20])
Top_25_in_otu <- intersect(Top_25, rownames(species_scores))
taxa_arrows <- as.data.frame(species_scores[Top_25_in_otu, ] * 1.8)
taxa_arrows$TaxonID <- rownames(taxa_arrows)

##Taxonomy labels##
tax_tbl <- as.data.frame(tax_table(physeq_raw))
tax_tbl_clean <- apply(tax_tbl, 2, function(x) gsub("^[a-zA-Z]__", "", x))
tax_tbl_clean <- as.data.frame(tax_tbl_clean)

taxa_arrows$Label <- sapply(taxa_arrows$TaxonID, function(tid) {
  family <- tax_tbl_clean[tid, "Family"]
  genus <- tax_tbl_clean[tid, "Genus"]
  order <- tax_tbl_clean[tid, "Order"]
  
  if (!is.na(family) && family != "" && !is.na(genus) && genus != "") {
    paste0(family, " ", genus)
  } else if (!is.na(family) && family != "") {
    paste0(family, " sp.")
  } else if (!is.na(order) && order != "") {
    paste0("Unclassified ", order)
  } else {
    "Unclassified"
  }
})
##ggplot2 colors and shapes##
myclos <- brewer.pal(n = 4, name = "RdBu")
shapes <- c(21, 22, 24, 25)

##Plot##
p <- ggplot(data = site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Cold_Storage, shape = Tissue), size = 4, color = "black") +
  scale_fill_manual(values = myclos) +
  scale_shape_manual(values = shapes[1:length(unique(site_scores$Tissue))]) +
  
##Top metabolite arrows##
  geom_segment(data = metabolite_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkblue", size = 1.2) +
  geom_text_repel(
    data = metabolite_arrows,
    aes(x = NMDS1, y = NMDS2, label = Metabolite),
    color = "darkblue", fontface = "bold", size = 4,
    max.overlaps = Inf) +
  
##Top taxa arrows##
geom_segment(data = taxa_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "gray30", size = 0.8) +
geom_text_repel(
    data = taxa_arrows,
    aes(x = NMDS1, y = NMDS2, label = Label),
    size = 3, fontface = "italic", color = "black",
    max.overlaps = Inf) +
  
theme_minimal() +
  labs(title = "Bacteria–Carbohydrates Acids NMDS Triplot",
       subtitle = paste("NMDS Stress =", round(ord.nmds.bray$stress, 3)),
       fill = "Cold Storage",
       shape = "Tissue Type") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

##Save plot##
ggsave("Bacteria_Carbohydrates_NMDS_Triplot.pdf", plot = p, width = 8, height = 7)


##Bacteria_KIWI vs Phenolics_Triplot##
##Load and scale metabolites##
Phenolics_Triplot <- read.table("Phenolics_Triplot.txt", header = TRUE, sep = "\t", row.names = 1)
Phenolics_Triplot_scaled <- scale(Phenolics_Triplot)

##Phyloseq object##
physeq_raw <- Bacteria_KIWI_2024_An_Merged_Glommed_100

##Build OTU matrix##
otu_mat <- as(otu_table(physeq_raw), "matrix")
if (taxa_are_rows(physeq_raw)) {
  otu_mat <- t(otu_mat)
}

##NMDS ordination##
ord.nmds.bray <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)

##Metadata##
meta <- as(sample_data(physeq_raw), "data.frame")
meta$SampleID <- rownames(meta)

##Extract NMDS site scores##
site_scores <- as.data.frame(scores(ord.nmds.bray, display = "sites"))
site_scores$SampleID <- rownames(site_scores)
site_scores <- left_join(site_scores, meta, by = "SampleID")

##Fit metabolites to ordination##
fit_env <- envfit(ord.nmds.bray, Phenolics_Triplot_scaled)

##Extract top metabolites based on R²##
r2_values <- fit_env$vectors$r
top5_idx <- order(r2_values, decreasing = TRUE)[1:6]
top5_arrows <- fit_env$vectors$arrows[top5_idx, , drop = FALSE]
metabolite_arrows <- as.data.frame(top5_arrows * 0.3)
metabolite_arrows$Metabolite <- rownames(top5_arrows)

##Extract taxa scores##
species_scores <- scores(ord.nmds.bray, display = "species")
Top_25 <- names(sort(taxa_sums(physeq_raw), decreasing = TRUE)[1:20])
Top_25_in_otu <- intersect(Top_25, rownames(species_scores))
taxa_arrows <- as.data.frame(species_scores[Top_25_in_otu, ] * 1.8)
taxa_arrows$TaxonID <- rownames(taxa_arrows)

##Taxonomy labels##
tax_tbl <- as.data.frame(tax_table(physeq_raw))
tax_tbl_clean <- apply(tax_tbl, 2, function(x) gsub("^[a-zA-Z]__", "", x))
tax_tbl_clean <- as.data.frame(tax_tbl_clean)

taxa_arrows$Label <- sapply(taxa_arrows$TaxonID, function(tid) {
  family <- tax_tbl_clean[tid, "Family"]
  genus <- tax_tbl_clean[tid, "Genus"]
  order <- tax_tbl_clean[tid, "Order"]
  
  if (!is.na(family) && family != "" && !is.na(genus) && genus != "") {
    paste0(family, " ", genus)
  } else if (!is.na(family) && family != "") {
    paste0(family, " sp.")
  } else if (!is.na(order) && order != "") {
    paste0("Unclassified ", order)
  } else {
    "Unclassified"
  }
})

##ggplot2 colors and shapes##
myclos <- brewer.pal(n = 4, name = "RdBu")
shapes <- c(21, 22, 24, 25)

##Plot##
p <- ggplot(data = site_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Cold_Storage, shape = Tissue), size = 4, color = "black") +
  scale_fill_manual(values = myclos) +
  scale_shape_manual(values = shapes[1:length(unique(site_scores$Tissue))]) +
  
## Top metabolite arrows##
geom_segment(data = metabolite_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkblue", size = 1.2) +
geom_text_repel(
    data = metabolite_arrows,
    aes(x = NMDS1, y = NMDS2, label = Metabolite),
    color = "darkblue", fontface = "bold", size = 4,
    max.overlaps = Inf) +
  
##Top taxa arrows##
geom_segment(data = taxa_arrows,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "gray30", size = 0.8) +
geom_text_repel(
    data = taxa_arrows,
    aes(x = NMDS1, y = NMDS2, label = Label),
    size = 3, fontface = "italic", color = "black",
    max.overlaps = Inf) +
  
theme_minimal() +
  labs(title = "Bacteria–Phenolics NMDS Triplot",
       subtitle = paste("NMDS Stress =", round(ord.nmds.bray$stress, 3)),
       fill = "Cold Storage",
       shape = "Tissue Type") +
theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 12))

##Save plot##
ggsave("Bacteria_Phenolics_NMDS_Triplot.pdf", plot = p, width = 8, height = 7)