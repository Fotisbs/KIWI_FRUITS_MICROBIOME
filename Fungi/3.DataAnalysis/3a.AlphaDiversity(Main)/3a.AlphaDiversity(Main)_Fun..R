#Alpha diversity Box Plots#
# Load required libraries
library(phyloseq)
library(reshape2)
library(agricolae)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Rarefy dataset
set.seed(123)
physeq_rare <- rarefy_even_depth(Fungi_KIWI_2024_Annotated, rngseed = 123, verbose = FALSE)

# Alpha diversity from rarefied data
adiv <- plot_richness(physeq_rare, 
                      measures = c("Observed", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))

# Combine alpha metrics
alpha_long <- adiv$data[, c("samples", "variable", "value")]

# Convert to wide format
alpha_wide <- dcast(alpha_long, samples ~ variable, value.var = "value")
row.names(alpha_wide) <- alpha_wide$samples

# Merge with metadata
alpha_wide_fact <- merge(alpha_wide, data.frame(sample_data(Fungi_KIWI_2024_Annotated)), 
                         by = "row.names")
row.names(alpha_wide_fact) <- alpha_wide_fact$Row.names
alpha_wide_fact <- alpha_wide_fact[ , !(colnames(alpha_wide_fact) %in% "Row.names")]

# Define alpha indices and grouping factors
mytestvars <- c("Observed", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
mytestfacts <- c("Tissue", "MCP", "Cold_Storage")

# Store stats
mystatsout <- list()

# Loop over grouping variables
for (mytestfact in mytestfacts) {
  cairo_pdf(paste0("alpha_div_boxplots_", mytestfact, "_smart_scaled.pdf"), height = 5, width = 9, onefile = TRUE)
  plot_list <- list()
  
  for (mytestvar in mytestvars) {
    
    myaovmatrix <- alpha_wide_fact
    myaovmatrix[[mytestfact]] <- factor(myaovmatrix[[mytestfact]])
    
    # Shapiro test
    shap_out <- by(myaovmatrix[[mytestvar]], myaovmatrix[[mytestfact]], shapiro.test)
    shap_pvals <- sapply(shap_out, function(x) x$p.value)
    mystatsout[[mytestvar]][[mytestfact]][["shapiro"]] <- shap_pvals
    
    non_normal <- any(shap_pvals < 0.05)
    test_type <- ifelse(non_normal, "Kruskal", "ANOVA")
    
    if (non_normal) {
      mykrusk <- kruskal(myaovmatrix[[mytestvar]], myaovmatrix[[mytestfact]], group = TRUE, p.adj = "BH")
      mystatsout[[mytestvar]][[mytestfact]][["kruskal"]] <- mykrusk
      group_table <- mykrusk$groups
    } else {
      myform <- as.formula(paste0("`", mytestvar, "` ~ ", mytestfact))
      mymod <- aov(myform, data = myaovmatrix)
      mysumaov <- summary(mymod)
      myHSDtest <- HSD.test(mymod, mytestfact, group = TRUE)
      mystatsout[[mytestvar]][[mytestfact]][["ANOVA"]] <- mysumaov
      mystatsout[[mytestvar]][[mytestfact]][["HSD test"]] <- myHSDtest
      group_table <- myHSDtest$groups
    }
    
    # Group label setup
    label_df <- data.frame(Group = rownames(group_table),
                           Letters = group_table$groups,
                           stringsAsFactors = FALSE)
    
    # Smart y-axis scaling
    range_y <- range(myaovmatrix[[mytestvar]], na.rm = TRUE)
    y_range <- diff(range_y)
    y_pad <- 0.1 * y_range
    y_limits <- c(range_y[1] - y_pad, range_y[2] + y_pad)
    
    if (mytestvar == "Simpson") {
      y_limits[1] <- max(0, y_limits[1])
      y_limits[2] <- min(1.05, y_limits[2])
    }
    
    custom_scale <- scale_y_continuous(limits = y_limits)
    label_df$y <- rep(y_limits[2], nrow(label_df))
    
    # Final plot
    p <- ggplot(myaovmatrix, aes(x = .data[[mytestfact]], 
                                 y = .data[[mytestvar]], 
                                 fill = .data[[mytestfact]])) +
      geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.8) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
      geom_text(data = label_df,
                aes(x = Group, y = y, label = Letters), 
                inherit.aes = FALSE, size = 4, fontface = "italic") +
      theme_minimal() +
      labs(title = paste0(mytestvar, " (", test_type, ")"),
           y = mytestvar, x = NULL) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 10)) +
      scale_fill_brewer(palette = "Set2") +
      custom_scale
    
    plot_list[[mytestvar]] <- p
  }
  
  grid.arrange(grobs = plot_list, ncol = 3)
  dev.off()
}