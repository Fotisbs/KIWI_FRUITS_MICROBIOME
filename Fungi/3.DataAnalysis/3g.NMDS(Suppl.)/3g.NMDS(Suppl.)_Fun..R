##Seperate Tissue Type Samples (Pericarp/Placenta)
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

###NMDS Pericarp Start###
Fungi_KIWI_2024_Pericarp_100
##ordination matrix
ord.nmds.bray1 <- ordinate(Fungi_KIWI_2024_Pericarp_100, method="NMDS", distance="bray")
#Cold Storage_MCP
plot_ordination(Fungi_KIWI_2024_Pericarp_100, ord.nmds.bray1, color="Cold_Storage", shape ="MCP", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
#MCP_Cold Storage
plot_ordination(Fungi_KIWI_2024_Pericarp_100, ord.nmds.bray1, color="MCP", shape ="Cold_Storage", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
###NMDS Pericarp END###

##Placenta
###NMDS Placenta Start##
Fungi_KIWI_2024_Placenta_100
##ordination matrix
ord.nmds.bray2 <- ordinate(Fungi_KIWI_2024_Placenta_100, method="NMDS", distance="bray")
#Cold Storage_MCP
plot_ordination(Fungi_KIWI_2024_Placenta_100, ord.nmds.bray2, color="Cold_Storage", shape ="MCP", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray2$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
#MCP_Cold Storage
plot_ordination(Fungi_KIWI_2024_Placenta_100, ord.nmds.bray2, color="MCP", shape ="Cold_Storage", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray2$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()
###NMDS Placenta END###