##NMDS Analysis#
##Relative abundance transformation 100%
Fungi_KIWI_2024_An_100 <- transform_sample_counts(Fungi_KIWI_2024_Annotated, function(OTU) 100*OTU/sum(OTU))

saveRDS(Fungi_KIWI_2024_An_100, file = "Fungi_KIWI_2024_An_100.RDS")

##NMDS Start##
##ordination matrix
ord.nmds.bray1 <- ordinate(Fungi_KIWI_2024_An_100, method="NMDS", distance="bray")

##One Factor Checked
#Tissue
plot_ordination(Fungi_KIWI_2024_An_100, ord.nmds.bray1, color="Tissue", shape ="", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00")) + theme_bw()

#MCP
plot_ordination(Fungi_KIWI_2024_An_100, ord.nmds.bray1, color="MCP", shape ="", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00")) + theme_bw()

#Cold Stor.
plot_ordination(Fungi_KIWI_2024_An_100, ord.nmds.bray1, color="Cold_Storage", shape ="", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C")) + theme_bw()

###NMDS END###