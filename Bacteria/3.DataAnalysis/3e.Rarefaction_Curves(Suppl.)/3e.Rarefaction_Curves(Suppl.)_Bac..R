##Rarefaction Curves##
Bacteria_KIWI_2024_Annotated

a <- data.frame((otu_table(Bacteria_KIWI_2024_Annotated)))

rarecurve(a, step=50, cex=0.5, label = F, xlim=c(0,1100), ylim =c(0,300))