##Fungi KIWI 2024 PERMANOVA##
##Relative abundance transformation 100%
Fungi_KIWI_2024_An_100 <- transform_sample_counts(Fungi_KIWI_2024_Annotated, function(OTU) 100*OTU/sum(OTU))

saveRDS(Fungi_KIWI_2024_An_100, file = "Fungi_KIWI_2024_An_100.RDS")

##PERMANOVA (Instead of CCA) Start##
mypermanova_Fungi_KIWI_2024_An_100 <- adonis2(Fungi_KIWI_2024_An_100@otu_table ~ MCP + Cold_Storage + Tissue, method = "bray", data = data.frame(Fungi_KIWI_2024_An_100@sam_data), by = "terms")
##PERMANOVA (Instead of CCA) END##

##pair-wise PERMANOVA Start##
##MCP
mycmpfactor_Fungi_KIWI_2024_An_100 <- interaction(data.frame(Fungi_KIWI_2024_An_100@sam_data)$MCP) 

mympairwiseperm_Fungi_KIWI_2024_An_100 <- pairwise.adonis(Fungi_KIWI_2024_An_100@otu_table, mycmpfactor_Fungi_KIWI_2024_An_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Cold Storage
mycmpfactor_Fungi_KIWI_2024_An_100 <- interaction(data.frame(Fungi_KIWI_2024_An_100@sam_data)$Cold_Storage) 

mympairwiseperm_Fungi_KIWI_2024_An_100 <- pairwise.adonis(Fungi_KIWI_2024_An_100@otu_table, mycmpfactor_Fungi_KIWI_2024_An_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Tissue
mycmpfactor_Fungi_KIWI_2024_An_100 <- interaction(data.frame(Fungi_KIWI_2024_An_100@sam_data)$Tissue) 

mympairwiseperm_Fungi_KIWI_2024_An_100 <- pairwise.adonis(Fungi_KIWI_2024_An_100@otu_table, mycmpfactor_Fungi_KIWI_2024_An_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
###pair-wise PERMANOVA END###

##Separate Tissue here (Pericarp/Placenta)
##Pericarp
Fungi_KIWI_2024_Pericarp <- subset_samples(Fungi_KIWI_2024_Annotated, Tissue=="Pericarp")

Fungi_KIWI_2024_Pericarp <- prune_taxa(taxa_sums(Fungi_KIWI_2024_Pericarp)>0,Fungi_KIWI_2024_Pericarp)

saveRDS(Fungi_KIWI_2024_Pericarp, file = "Fungi_KIWI_2024_Pericarp.RDS")

Fungi_KIWI_2024_Pericarp

##Placenta
Fungi_KIWI_2024_Placenta <- subset_samples(Fungi_KIWI_2024_Annotated, Tissue=="Placenta")

Fungi_KIWI_2024_Placenta <- prune_taxa(taxa_sums(Fungi_KIWI_2024_Placenta)>0,Fungi_KIWI_2024_Placenta)

saveRDS(Fungi_KIWI_2024_Placenta, file = "Fungi_KIWI_2024_Placenta.RDS")

Fungi_KIWI_2024_Placenta

##RA % per Tissue
Fungi_KIWI_2024_Pericarp
Fungi_KIWI_2024_Placenta

Fungi_KIWI_2024_Pericarp_100 <- transform_sample_counts(Fungi_KIWI_2024_Pericarp, function(OTU) 100*OTU/sum(OTU))

saveRDS(Fungi_KIWI_2024_Pericarp_100, file = "Fungi_KIWI_2024_Pericarp_100.RDS")

Fungi_KIWI_2024_Placenta_100 <- transform_sample_counts(Fungi_KIWI_2024_Placenta, function(OTU) 100*OTU/sum(OTU))

saveRDS(Fungi_KIWI_2024_Placenta_100, file = "Fungi_KIWI_2024_Placenta_100.RDS")

Fungi_KIWI_2024_Pericarp_100
Fungi_KIWI_2024_Placenta_100

##PERMANOVA Pericarp Start ##
Fungi_KIWI_2024_Pericarp_100

mypermanova_Fungi_KIWI_2024_Pericarp_100 <- adonis2(Fungi_KIWI_2024_Pericarp_100@otu_table ~ MCP + Cold_Storage, method = "bray", data = data.frame(Fungi_KIWI_2024_Pericarp_100@sam_data), by = "terms")
##PERMANOVA Pericarp  END##

###pair-wise PERMANOVA Pericarp Start###
#MCP
mycmpfactor_Fungi_KIWI_2024_Pericarp_100 <- interaction(data.frame(Fungi_KIWI_2024_Pericarp_100@sam_data)$MCP) 

mympairwiseperm_Fungi_KIWI_2024_Pericarp_100 <- pairwise.adonis(Fungi_KIWI_2024_Pericarp_100@otu_table, mycmpfactor_Fungi_KIWI_2024_Pericarp_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

#Cold Storage
mycmpfactor_Fungi_KIWI_2024_Pericarp_100 <- interaction(data.frame(Fungi_KIWI_2024_Pericarp_100@sam_data)$Cold_Storage) 

mympairwiseperm_Fungi_KIWI_2024_Pericarp_100 <- pairwise.adonis(Fungi_KIWI_2024_Pericarp_100@otu_table, mycmpfactor_Fungi_KIWI_2024_Pericarp_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
###pair-wise PERMANOVA Pericarp END###

##Placenta
##PERMANOVA Placenta Start ##
Fungi_KIWI_2024_Placenta_100

mypermanova_Fungi_KIWI_2024_Placenta_100 <- adonis2(Fungi_KIWI_2024_Placenta_100@otu_table ~ Cold_Storage + MCP, method = "bray", data = data.frame(Fungi_KIWI_2024_Placenta_100@sam_data), by = "terms")
##PERMANOVA Placenta  END##

###pair-wise PERMANOVA Placenta Start###
#MCP
mycmpfactor_Fungi_KIWI_2024_Placenta_100 <- interaction(data.frame(Fungi_KIWI_2024_Placenta_100@sam_data)$MCP) 

mympairwiseperm_Fungi_KIWI_2024_Placenta_100 <- pairwise.adonis(Fungi_KIWI_2024_Placenta_100@otu_table, mycmpfactor_Fungi_KIWI_2024_Placenta_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

#Cold Storage
mycmpfactor_Fungi_KIWI_2024_Placenta_100 <- interaction(data.frame(Fungi_KIWI_2024_Placenta_100@sam_data)$Cold_Storage)

mympairwiseperm_Fungi_KIWI_2024_Placenta_100 <- pairwise.adonis(Fungi_KIWI_2024_Placenta_100@otu_table, mycmpfactor_Fungi_KIWI_2024_Placenta_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
###pair-wise PERMANOVA Pericarp END###