# ============================================================================
# extract environmental factors
# ============================================================================
library(geodata)
library(raster)
library(sp)
library(rgdal)

# get soil grid layers (tif files)
gph <- soil_world(var="phh2o", depth=5, path="soilgrid")
gni <- soil_world(var="nitrogen", depth=5, path="soilgrid")
gsoc <- soil_world(var="soc", depth=5, path="soilgrid")
gsand <- soil_world(var="sand", depth=5, path="soilgrid")
gsilt <- soil_world(var="silt", depth=5, path="soilgrid")
gclay <- soil_world(var="clay", depth=5, path="soilgrid")
ph <- raster(gph)
ni <- raster(gni)
soc <- raster(gsoc)
sand <- raster(gsand)
silt <- raster(gsilt)
clay <- raster(gclay)

# get WorldClim and evapo layers,
# download from https://www.worldclim.org/data/worldclim21.html and
# https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/2
data_path1<-dir("database/worldclim",pattern = "tif",full.names=TRUE)
worldclim=stack(data_path)
data_path2<-dir("database/Global Aridity and PET Database",pattern = "tif",full.names=TRUE)
evapo=raster(data_path2[1])

#change the filename as your own location file
Samples<-read.table("location.txt",row.names = 1,header =T ,sep = "\t") 
#choose worldclim as the example
EF_data<-extract(worldclim,Samples,method="bilinear")

# ============================================================================
# CONCISE DBRDA ANALYSIS WITH AUTO-EXTRACTED STATISTICS
# ============================================================================

# Load libraries
library(ggplot2)
library(patchwork)
library(vegan)
library(missForest)
library(dplyr)
library(tibble)
library(ggrepel)
library(showtext)

# Setup
font_add_google("Tinos", "myFont")
showtext_auto()

# ============================================================================
# DATA LOADING
# ============================================================================

# Load data
files <- c("microcoleus_env0213.txt", "og_dist.matrix", "pro_dist.matrix", 
           "gene_dist.matrix", "ani_dist.matrix")
names <- c("Mv.env", "og_dist", "pro_dist", "gene_dist", "ani_dist")

for(i in seq_along(files)) {
  assign(names[i], read.table(files[i], header=TRUE, row.names=1, sep="\t", check.names=FALSE))
}

# Impute missing values
char_cols <- sapply(Mv.env, is.character)
if(any(char_cols)) Mv.env[char_cols] <- lapply(Mv.env[char_cols], as.factor)
Mv.env.imputed <- missForest(Mv.env)$ximp

# ============================================================================
# ANALYSIS CONFIGURATION
# ============================================================================

# Check available columns
cat("Available columns:", paste(colnames(Mv.env.imputed), collapse=", "), "\n")

# Configuration - adjust env_vars based on your actual column names
config <- list(
  ani = list(dist=ani_dist, env_vars=c("precipitation", "evapo", "AI", "WHC", "pH", "radiation"), 
             title="Whole genome (ANI)", pos=c(-1, 2.5)),
  og = list(dist=og_dist, env_vars=c("precipitation", "evapo", "AI", "WHC", "pH", "radiation"), 
            title="Orthologous genes", pos=c(-0.7, 0.9)),
  gene = list(dist=gene_dist, env_vars=c("precipitation", "evapo", "AI", "WHC", "pH", "radiation", "elevation"), 
              title="Core genes", pos=c(-0.7, 0.6)),
  pro = list(dist=pro_dist, env_vars=c("precipitation", "evapo", "AI", "WHC", "pH", "radiation", "elevation"), 
             title="Core proteins", pos=c(-0.7, 0.6))
)

# ============================================================================
# ANALYSIS FUNCTION
# ============================================================================

run_dbrda <- function(cfg, name) {
  cat("\nAnalyzing:", name, "\n")
  
  # Check if columns exist
  missing <- setdiff(cfg$env_vars, colnames(Mv.env.imputed))
  if(length(missing) > 0) {
    cat("Missing columns:", paste(missing, collapse=", "), "\n")
    return(NULL)
  }
  
  # Subset data
  env_data <- Mv.env.imputed[, cfg$env_vars, drop=FALSE]
  
  # Run dbRDA
  model <- dbrda(cfg$dist ~ ., data=env_data)
  
  # Extract statistics
  anova_res <- anova(model)
  adonis_res <- adonis2(cfg$dist ~ ., data=env_data, permutations=999)
  
  # Calculate variance explained
  eigenvals <- model$CCA$eig
  total_var <- sum(eigenvals)
  axis1_pct <- round((eigenvals[1]/total_var)*100, 1)
  axis2_pct <- round((eigenvals[2]/total_var)*100, 1)
  
  # Extract scores
  site_scores <- as.data.frame(scores(model, display="sites")) %>%
    merge(Mv.env.imputed["phylogroup"], by="row.names") %>%
    column_to_rownames("Row.names")
  
  biplot_data <- as.data.frame(scores(model, display="bp")) %>%
    rownames_to_column("variable")
  
  # Calculate hulls
  hulls <- site_scores %>%
    group_by(phylogroup) %>%
    slice(chull(dbRDA1, dbRDA2)) %>%
    ungroup()
  
  # Compile statistics
  stats <- list(
    F = round(anova_res$F[1], 2),
    p = round(anova_res$`Pr(>F)`[1], 3),
    R2 = round(adonis_res$R2[1], 3),
    axis1 = axis1_pct,
    axis2 = axis2_pct
  )
  
  cat(sprintf("  F=%.2f, p=%.3f, R²=%.3f, Axis1=%.1f%%, Axis2=%.1f%%\n",
              stats$F, stats$p, stats$R2, stats$axis1, stats$axis2))
  
  return(list(site_scores=site_scores, biplot_data=biplot_data, 
              hulls=hulls, stats=stats))
}

# ============================================================================
# PLOTTING FUNCTION
# ============================================================================

create_plot <- function(results, cfg, show_legend=FALSE) {
  if(is.null(results)) return(NULL)
  
  ggplot(results$site_scores, aes(dbRDA1, dbRDA2)) +
    # Hulls
    geom_polygon(data=results$hulls, aes(group=phylogroup), 
                 alpha=0.3, fill="grey80", color="grey60", linewidth=0.1, show.legend=FALSE) +
    
    # Points
    geom_point(aes(fill=phylogroup), shape=21, size=2.5, color="white", 
               stroke=0.01, alpha=0.7) +
    
    # Arrows
    geom_segment(data=results$biplot_data, aes(x=0, y=0, xend=dbRDA1, yend=dbRDA2),
                 arrow=arrow(length=unit(0.12,"cm")), color="navy", linewidth=0.4) +
    
    # Labels
    geom_text_repel(data=results$biplot_data, aes(x=dbRDA1, y=dbRDA2, label=variable),
                    color="navy", size=2, family="myFont", max.overlaps=Inf) +
    
    # Statistics annotation
    annotate("text", x=cfg$pos[1], y=cfg$pos[2],
             label=sprintf("atop(paste('F = %.2f, p = %.3f'), paste('R'^2*' = %.3f'))",
                           results$stats$F, results$stats$p, results$stats$R2),
             parse=TRUE, size=1.8, family="myFont") +
    
    # Reference lines
    geom_vline(xintercept=0, color="dimgrey", linewidth=0.2, linetype=2) +
    geom_hline(yintercept=0, color="dimgrey", linewidth=0.2, linetype=2) +
    
    # Scales
    scale_fill_manual(values=c("#756BB1", "#31A354", "#BCBDDC", "#3182BD", 
                               "#FDAE6B", "#A1D99B", "#9ECAE1")) +
    scale_x_continuous(name=sprintf("dbRDA 1 (%.1f%%)", results$stats$axis1)) +
    scale_y_continuous(name=sprintf("dbRDA 2 (%.1f%%)", results$stats$axis2)) +
    
    # Theme
    ggtitle(cfg$title) +
    theme_bw() +
    theme(
      axis.title=element_text(size=7, family="myFont"),
      axis.text=element_text(size=5, family="myFont"),
      panel.grid=element_blank(),
      plot.title=element_text(size=8, hjust=0.5, family="myFont", face="bold"),
      legend.position=if(show_legend) "right" else "none",
      legend.text=element_text(size=6, family="myFont"),
      legend.title=element_blank(),
      legend.key.size=unit(0.4,"cm")
    )
}

# ============================================================================
# RUN ANALYSES
# ============================================================================

results <- list()
plots <- list()

for(name in names(config)) {
  results[[name]] <- run_dbrda(config[[name]], name)
  if(!is.null(results[[name]])) {
    plots[[name]] <- create_plot(results[[name]], config[[name]], show_legend=TRUE)
  }
}

# ============================================================================
# COMBINE AND SAVE
# ============================================================================

if(length(plots) > 0) {
  # Create summary table
  summary_df <- data.frame(
    Analysis = names(results)[!sapply(results, is.null)],
    F_statistic = sapply(results[!sapply(results, is.null)], function(x) x$stats$F),
    p_value = sapply(results[!sapply(results, is.null)], function(x) x$stats$p),
    R_squared = sapply(results[!sapply(results, is.null)], function(x) x$stats$R2),
    dbRDA1_percent = sapply(results[!sapply(results, is.null)], function(x) x$stats$axis1),
    dbRDA2_percent = sapply(results[!sapply(results, is.null)], function(x) x$stats$axis2)
  )
  
  print("ANALYSIS SUMMARY:")
  print(summary_df)
  write.csv(summary_df, "dbrda_summary.csv", row.names=FALSE)
  
  # Combine plots
  if(length(plots) == 4) {
    final_plot <- wrap_plots(plots, ncol=2, guides="collect") &
      theme(legend.position="right", legend.text=element_text(size=8, family="myFont"))
  } else {
    final_plot <- wrap_plots(plots, guides="collect") &
      theme(legend.position="right", legend.text=element_text(size=8, family="myFont"))
  }
  
  # Add title
  final_plot <- final_plot + 
    plot_annotation(title="dbRDA Analysis - Comparative Genomics",
                    theme=theme(plot.title=element_text(size=12, hjust=0.5, family="myFont")))
  
  # Display and save
  print(final_plot)
  ggsave("dbrda_analysis.pdf", final_plot, width=18, height=16, units="cm", device=cairo_pdf)
  
  cat("\nFiles saved:\n")
  cat("- Combined plot: dbrda_analysis.pdf\n")
  cat("- Summary table: dbrda_summary.csv\n")
  
} else {
  cat("No successful analyses. Check your column names.\n")
}

# Cleanup
showtext_auto(FALSE)
cat("Analysis completed!\n")
