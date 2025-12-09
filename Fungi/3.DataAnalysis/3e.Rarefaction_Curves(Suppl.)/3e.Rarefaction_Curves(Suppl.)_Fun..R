##Rarefaction Curves##
Fungi_KIWI_2024_Annotated

a <- data.frame((otu_table(Fungi_KIWI_2024_Annotated)))

rarecurve(a, step=50, cex=0.5, label = F, xlim=c(0,1100), ylim =c(0,100))