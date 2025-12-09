##Bacteria KIWI 2024 PERMANOVA##
##Relative abundance transformation 100%
Bacteria_KIWI_2024_An_100 <- transform_sample_counts(Bacteria_KIWI_2024_Annotated, function(OTU) 100*OTU/sum(OTU))

saveRDS(Bacteria_KIWI_2024_An_100, file = "Bacteria_KIWI_2024_An_100.RDS")

##PERMANOVA (Instead of CCA) Start##
mypermanova_Bacteria_KIWI_2024_An_100 <- adonis2(Bacteria_KIWI_2024_An_100@otu_table ~ MCP + Cold_Storage + Tissue, method = "bray", data = data.frame(Bacteria_KIWI_2024_An_100@sam_data), by = "terms")
##PERMANOVA (Instead of CCA) END##

##pair-wise PERMANOVA Start##
##MCP
mycmpfactor_Bacteria_KIWI_2024_An_100 <- interaction(data.frame(Bacteria_KIWI_2024_An_100@sam_data)$MCP) 

mympairwiseperm_Bacteria_KIWI_2024_An_100 <- pairwise.adonis(Bacteria_KIWI_2024_An_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_An_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Cold Storage
mycmpfactor_Bacteria_KIWI_2024_An_100 <- interaction(data.frame(Bacteria_KIWI_2024_An_100@sam_data)$Cold_Storage) 

mympairwiseperm_Bacteria_KIWI_2024_An_100 <- pairwise.adonis(Bacteria_KIWI_2024_An_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_An_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Tissue
mycmpfactor_Bacteria_KIWI_2024_An_100 <- interaction(data.frame(Bacteria_KIWI_2024_An_100@sam_data)$Tissue) 

mympairwiseperm_Bacteria_KIWI_2024_An_100 <- pairwise.adonis(Bacteria_KIWI_2024_An_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_An_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
###pair-wise PERMANOVA END###

##Separate Tissue here (Pericarp/Placenta)
##Pericarp
Bacteria_KIWI_2024_Pericarp <- subset_samples(Bacteria_KIWI_2024_Annotated, Tissue=="Pericarp")

Bacteria_KIWI_2024_Pericarp <- prune_taxa(taxa_sums(Bacteria_KIWI_2024_Pericarp)>0,Bacteria_KIWI_2024_Pericarp)

saveRDS(Bacteria_KIWI_2024_Pericarp, file = "Bacteria_KIWI_2024_Pericarp.RDS")

Bacteria_KIWI_2024_Pericarp

##Placenta
Bacteria_KIWI_2024_Placenta <- subset_samples(Bacteria_KIWI_2024_Annotated, Tissue=="Placenta")

Bacteria_KIWI_2024_Placenta <- prune_taxa(taxa_sums(Bacteria_KIWI_2024_Placenta)>0,Bacteria_KIWI_2024_Placenta)

saveRDS(Bacteria_KIWI_2024_Placenta, file = "Bacteria_KIWI_2024_Placenta.RDS")

Bacteria_KIWI_2024_Placenta

##RA % per Tissue
Bacteria_KIWI_2024_Pericarp
Bacteria_KIWI_2024_Placenta

Bacteria_KIWI_2024_Pericarp_100 <- transform_sample_counts(Bacteria_KIWI_2024_Pericarp, function(OTU) 100*OTU/sum(OTU))

saveRDS(Bacteria_KIWI_2024_Pericarp_100, file = "Bacteria_KIWI_2024_Pericarp_100.RDS")

Bacteria_KIWI_2024_Placenta_100 <- transform_sample_counts(Bacteria_KIWI_2024_Placenta, function(OTU) 100*OTU/sum(OTU))

saveRDS(Bacteria_KIWI_2024_Placenta_100, file = "Bacteria_KIWI_2024_Placenta_100.RDS")

Bacteria_KIWI_2024_Pericarp_100
Bacteria_KIWI_2024_Placenta_100

##PERMANOVA Pericarp Start ##
Bacteria_KIWI_2024_Pericarp_100

mypermanova_Bacteria_KIWI_2024_Pericarp_100 <- adonis2(Bacteria_KIWI_2024_Pericarp_100@otu_table ~ MCP + Cold_Storage, method = "bray", data = data.frame(Bacteria_KIWI_2024_Pericarp_100@sam_data), by = "terms")
##PERMANOVA Pericarp  END##

###pair-wise PERMANOVA Pericarp Start###
#MCP
mycmpfactor_Bacteria_KIWI_2024_Pericarp_100 <- interaction(data.frame(Bacteria_KIWI_2024_Pericarp_100@sam_data)$MCP) 

mympairwiseperm_Bacteria_KIWI_2024_Pericarp_100 <- pairwise.adonis(Bacteria_KIWI_2024_Pericarp_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_Pericarp_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

#Cold Storage
mycmpfactor_Bacteria_KIWI_2024_Pericarp_100 <- interaction(data.frame(Bacteria_KIWI_2024_Pericarp_100@sam_data)$Cold_Storage) 

mympairwiseperm_Bacteria_KIWI_2024_Pericarp_100 <- pairwise.adonis(Bacteria_KIWI_2024_Pericarp_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_Pericarp_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
###pair-wise PERMANOVA Pericarp END###

##Placenta
##PERMANOVA Placenta Start ##
Bacteria_KIWI_2024_Placenta_100

mypermanova_Bacteria_KIWI_2024_Placenta_100 <- adonis2(Bacteria_KIWI_2024_Placenta_100@otu_table ~ Cold_Storage + MCP, method = "bray", data = data.frame(Bacteria_KIWI_2024_Placenta_100@sam_data), by = "terms")
##PERMANOVA Placenta  END##

###pair-wise PERMANOVA Placenta Start###
#MCP
mycmpfactor_Bacteria_KIWI_2024_Placenta_100 <- interaction(data.frame(Bacteria_KIWI_2024_Placenta_100@sam_data)$MCP) 

mympairwiseperm_Bacteria_KIWI_2024_Placenta_100 <- pairwise.adonis(Bacteria_KIWI_2024_Placenta_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_Placenta_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

#Cold Storage
mycmpfactor_Bacteria_KIWI_2024_Placenta_100 <- interaction(data.frame(Bacteria_KIWI_2024_Placenta_100@sam_data)$Cold_Storage)

mympairwiseperm_Bacteria_KIWI_2024_Placenta_100 <- pairwise.adonis(Bacteria_KIWI_2024_Placenta_100@otu_table, mycmpfactor_Bacteria_KIWI_2024_Placenta_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
###pair-wise PERMANOVA Pericarp END###