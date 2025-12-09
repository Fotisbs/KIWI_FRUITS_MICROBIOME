###Bar Plots Start###
##top 20 taxa at Genus level##
Fungi_KIWI_2024_Annotated

Fungi_KIWI_2024_An_100 <- transform_sample_counts(Fungi_KIWI_2024_Annotated, function(OTU) 100*OTU/sum(OTU))

Fungi_KIWI_2024_100_G.G <- tax_glom(Fungi_KIWI_2024_An_100, taxrank = "Genus")

Fungi_KIWI_2024_100_G.G

###Genus Start###
myTaxa11_Fungi_KIWI_2024_100 <- names(sort(taxa_sums(Fungi_KIWI_2024_100_G.G), decreasing = TRUE)[1:20])  

Top11_Fungi_KIWI_2024_100 <- prune_taxa(myTaxa11_Fungi_KIWI_2024_100, Fungi_KIWI_2024_100_G.G)

taxa_names(Top11_Fungi_KIWI_2024_100)

mytax11_Fungi_KIWI_2024_100 <- data.frame(tax_table(Top11_Fungi_KIWI_2024_100), stringsAsFactors = F)

# For ITS - Remove letter from taxonomy
for (i in c(1:nrow(mytax11_Fungi_KIWI_2024_100))) {
  for(j in c(1:ncol(mytax11_Fungi_KIWI_2024_100))) {
    mytax11_Fungi_KIWI_2024_100[i,j] <- gsub("[a-z]__","",mytax11_Fungi_KIWI_2024_100[i,j])
  }
}

mytxplot11_Fungi_KIWI_2024_100 <- data.frame(OTU = row.names(mytax11_Fungi_KIWI_2024_100), 
                                             txplt = paste(row.names(mytax11_Fungi_KIWI_2024_100), " ", mytax11_Fungi_KIWI_2024_100$Genus,  ":", mytax11_Fungi_KIWI_2024_100$Genus,  sep = ""))

row.names(mytxplot11_Fungi_KIWI_2024_100) <- mytxplot11_Fungi_KIWI_2024_100$OTU

taxa_names(Top11_Fungi_KIWI_2024_100) <- mytxplot11_Fungi_KIWI_2024_100[taxa_names(Top11_Fungi_KIWI_2024_100),"txplt"]

pdf(file = "L.pdf", width = 14, height = 7)

mycols <- c("burlywood2","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#E41A1C","#377EB8","coral2","#984EA3","#FF7F00","#A6CEE3","purple","#B2DF8A","aquamarine3","#FB9A99","blue2","red4")


plot_bar(Top11_Fungi_KIWI_2024_100, x="Cold_Storage_Rep", fill="Genus", title = "") + facet_grid(cols = vars(Tissue),rows = vars(MCP)) + geom_col() + scale_fill_manual(values = mycols)
###Genus END###