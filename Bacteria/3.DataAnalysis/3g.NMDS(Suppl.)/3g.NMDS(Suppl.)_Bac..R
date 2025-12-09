##Seperate Tissue Type Samples (Pericarp/Placenta)
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

###NMDS Pericarp Start###
Bacteria_KIWI_2024_Pericarp_100
##ordination matrix
ord.nmds.bray1 <- ordinate(Bacteria_KIWI_2024_Pericarp_100, method="NMDS", distance="bray")
#Cold Storage_MCP
plot_ordination(Bacteria_KIWI_2024_Pericarp_100, ord.nmds.bray1, color="Cold_Storage", shape ="MCP", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
#MCP_Cold Storage
plot_ordination(Bacteria_KIWI_2024_Pericarp_100, ord.nmds.bray1, color="MCP", shape ="Cold_Storage", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
###NMDS Pericarp END###

##Placenta
###NMDS Placenta Start##
Bacteria_KIWI_2024_Placenta_100
##ordination matrix
ord.nmds.bray2 <- ordinate(Bacteria_KIWI_2024_Placenta_100, method="NMDS", distance="bray")
#Cold Storage_MCP
plot_ordination(Bacteria_KIWI_2024_Placenta_100, ord.nmds.bray2, color="Cold_Storage", shape ="MCP", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray2$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
#MCP_Cold Storage
plot_ordination(Bacteria_KIWI_2024_Placenta_100, ord.nmds.bray2, color="MCP", shape ="Cold_Storage", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray2$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
###NMDS Placenta END###