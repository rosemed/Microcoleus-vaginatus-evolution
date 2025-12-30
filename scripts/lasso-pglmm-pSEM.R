#install.packages(c("ape", "phyr", "brms", "tidybayes", "ggplot2"))
library(ape)
library(phyr)
library(brms)
library(tidybayes)
library(ggplot2)
library(glmnet)
#library(glmmLasso)
library(lme4)
library(missForest)
library(tidyverse)
library(piecewiseSEM)

# Load and prepare data
data <- read.csv("SEM_data.csv", stringsAsFactors = FALSE)
cat("Dataset:", nrow(data), "observations ×", ncol(data), "variables\n")

# Convert non-numeric variables to factors
numeric_vars <- sapply(data, is.numeric)
data[!numeric_vars] <- lapply(data[!numeric_vars], as.factor)

# Handle missing data with missForest
set.seed(123)
data_imputed <- missForest(data, verbose = FALSE)$ximp
cat("✓ Missing data imputed\n")

# Standardize all numeric variables
data_standardized <- data_imputed
data_standardized[numeric_vars] <- scale(data_imputed[numeric_vars])

# Define variable sets
response_vars <- c("sel.gain", "SUR", "Ne")
abiotic_factors <- c("MAT", "precipitation", "elevation",
                      "AI", "WHC", "N", "pH", "SOC")

# Prepare complete dataset
all_vars <- c(response_vars, abiotic_factors)
available_vars <- intersect(all_vars, colnames(data_standardized))
complete_data <- data_standardized[complete.cases(data_standardized[available_vars]), ]
abiotic_available <- intersect(abiotic_factors, colnames(complete_data))


tree <- ape::read.tree("132_singlecopygene_reroot.treefile")
tree <- keep.tip(tree, intersect(tree$tip.label, data$strain))

dist=read.table("ani_dist.matrix.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
dist_mat=dist[rownames(dist) %in% complete_data$strain,colnames(dist) %in% complete_data$strain]

group_df=complete_data[,c(1,3)]
# Convert distance matrix to long format
dist_long <- as.data.frame(dist_mat) %>%
  rownames_to_column("strain1") %>%
  pivot_longer(-strain1, names_to = "strain2", values_to = "distance") %>%
  # Remove self-comparisons (where strain1 equals strain2)
  filter(strain1 != strain2)

# Merge with group information
dist_with_group <- dist_long %>%
  # Add group information for first strain
  left_join(group_df, by = c("strain1" = "strain")) %>%
  rename(phylogroup1 = phylogroup) %>%
  # Add group information for second strain
  left_join(group_df, by = c("strain2" = "strain")) %>%
  rename(phylogroup2 = phylogroup) %>%
  # Classify comparison type (within same group or between different groups)
  mutate(comparison = ifelse(phylogroup1 == phylogroup2, "within", "between"))

# Calculate average within-group distances for each phylogenetic group
within_group <- dist_with_group %>%
  filter(comparison == "within") %>%
  group_by(phylogroup1) %>%
  summarise(mean_within = mean(distance), .groups = "drop") %>%
  rename(phylogroup = phylogroup1)

# Calculate average between-group distances for each phylogenetic group
between_group <- dist_with_group %>%
  filter(comparison == "between") %>%
  # Collect comparisons for each group from both perspectives
  pivot_longer(cols = c(phylogroup1, phylogroup2), 
               names_to = "group_role", 
               values_to = "phylogroup") %>%
  group_by(phylogroup) %>%
  summarise(mean_between = mean(distance), .groups = "drop")

# Calculate aggregation index (within-group mean / between-group mean)
group_agg_index <- within_group %>%
  full_join(between_group, by = "phylogroup") %>%
  mutate(aggregation_index = mean_within / mean_between)

complete_data=complete_data %>% left_join(group_agg_index[,c(1,4)],by="phylogroup")

# ===============================================================================
# ENHANCED LASSO VARIABLE SELECTION
# ===============================================================================

cat("\n=== ENHANCED LASSO VARIABLE SELECTION ===\n")

# Create output directory
if (!dir.exists("lasso_outputs")) dir.create("lasso_outputs")

# Enhanced LASSO function with detailed diagnostics
perform_lasso_detailed <- function(response, predictors, data) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("DETAILED LASSO ANALYSIS FOR:", toupper(response), "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Prepare matrices
  X <- as.matrix(data[predictors])
  y <- data[[response]]
  
  cat("→ Dataset: n =", nrow(X), ", p =", ncol(X), "\n")
  cat("→ Predictors:", paste(predictors, collapse = ", "), "\n")
  cat("→ Response range:", round(min(y), 3), "to", round(max(y), 3), "\n\n")
  
  # Cross-validated LASSO
  set.seed(123)
  cat("→ Running 10-fold cross-validation...\n")
  cv_lasso <- cv.glmnet(X, y, alpha = 1, nfolds = 10, standardize = FALSE)
  
  # Extract coefficients
  coef_min <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  coef_1se <- as.matrix(coef(cv_lasso, s = "lambda.1se"))
  
  # Conservative selection (λ.1se)
  selected_1se <- rownames(coef_1se)[coef_1se != 0][-1]
  selected_min <- rownames(coef_min)[coef_min != 0][-1]
  
  # Detailed diagnostics
  cat("\n--- LASSO REGULARIZATION ANALYSIS ---\n")
  cat("λ.min  =", sprintf("%.6f", cv_lasso$lambda.min), 
      "(CV error =", sprintf("%.4f", min(cv_lasso$cvm)), ")\n")
  cat("λ.1se  =", sprintf("%.6f", cv_lasso$lambda.1se), 
      "(CV error =", sprintf("%.4f", cv_lasso$cvm[cv_lasso$lambda == cv_lasso$lambda.1se]), ")\n")
  
  cat("\n--- VARIABLE SELECTION COMPARISON ---\n")
  cat("Variables at λ.min :", length(selected_min), 
      ifelse(length(selected_min) > 0, paste("(", paste(selected_min, collapse = ", "), ")", sep = ""), ""), "\n")
  cat("Variables at λ.1se :", length(selected_1se), 
      ifelse(length(selected_1se) > 0, paste("(", paste(selected_1se, collapse = ", "), ")", sep = ""), ""), "\n")
  
  # Create coefficient table
  all_vars <- rownames(coef_1se)
  coef_detailed <- data.frame(
    Variable = all_vars,
    Coef_lambdamin = as.numeric(coef_min[all_vars, 1]),
    Coef_lambda1se = as.numeric(coef_1se[all_vars, 1]),
    Selected_lambdamin = all_vars %in% c("(Intercept)", selected_min),
    Selected_lambda1se = all_vars %in% c("(Intercept)", selected_1se),
    Abs_coef_1se = abs(as.numeric(coef_1se[all_vars, 1])),
    Correlation = NA,
    Correlation_SE = NA,
    stringsAsFactors = FALSE
  )
  
  # Add correlations and their standard errors
  for (i in 1:nrow(coef_detailed)) {
    var_name <- coef_detailed$Variable[i]
    if (var_name != "(Intercept)" && var_name %in% colnames(data)) {
      # Calculate correlation
      r_value <- cor(data[[var_name]], data[[response]], use = "complete.obs")
      coef_detailed$Correlation[i] <- r_value
      
      # Calculate standard error of correlation: SE(r) = sqrt((1-r²)/(n-2))
      if (!is.na(r_value)) {
        # Get sample size for complete observations
        complete_obs <- complete.cases(data[[var_name]], data[[response]])
        n <- sum(complete_obs)
        
        if (n > 2) {
          se_r <- sqrt((1 - r_value^2) / (n - 2))
          coef_detailed$Correlation_SE[i] <- se_r
        }
      }
    }
  }
  
  # Sort by importance
  coef_detailed <- coef_detailed[order(-coef_detailed$Abs_coef_1se), ]
  
  # Display top variables
  cat("\n--- TOP 10 VARIABLES BY |COEFFICIENT| AT λ.1se ---\n")
  display_vars <- head(coef_detailed[coef_detailed$Variable != "(Intercept)", ], 10)
  for (i in 1:min(10, nrow(display_vars))) {
    var_info <- display_vars[i, ]
    status_1se <- ifelse(var_info$Selected_lambda1se, "SELECTED", "excluded")
    status_min <- ifelse(var_info$Selected_lambdamin, "would select", "would exclude")
    
    # Format correlation with standard error
    if (!is.na(var_info$Correlation_SE)) {
      cor_text <- sprintf("r=%6.3f(±%5.3f)", var_info$Correlation, var_info$Correlation_SE)
    } else {
      cor_text <- sprintf("r=%6.3f", var_info$Correlation)
    }
    
    cat(sprintf("%2d. %-12s: λ.1se coef=%8.4f (%s), λ.min coef=%8.4f (%s), %s\n",
                i, var_info$Variable, var_info$Coef_lambda1se, status_1se,
                var_info$Coef_lambdamin, status_min, cor_text))
  }
  
  # Export results
  write.csv(coef_detailed, paste0("lasso_outputs/", response, "_detailed_coefficients.csv"), row.names = FALSE)
  
  cat("\n--- SELECTION SUMMARY ---\n")
  cat("CONSERVATIVE SELECTION (λ.1se):", length(selected_1se), "variables\n")
  if (length(selected_1se) > 0) {
    cat("Selected:", paste(selected_1se, collapse = ", "), "\n")
  } else {
    cat("No variables selected (intercept-only)\n")
  }
  
  return(list(selected = selected_1se, cv_model = cv_lasso, coefficients = coef_detailed))
}
# ===============================================================================
# HIERARCHICAL LASSO SELECTION
# ===============================================================================

lasso_results <- list()

# Step 1: sel.gain
cat("\n>>> STEP 1: sel.gain <<<")
lasso_results$sel.gain <- perform_lasso_detailed("sel.gain", c(abiotic_available,"aggregation_index"), complete_data)

# Step 2: SUR
cat("\n>>> STEP 2: SUR <<<")
lasso_results$SUR <- perform_lasso_detailed("SUR", c(abiotic_available,"aggregation_index"), complete_data)

# Step 3: Ne
cat("\n>>> STEP 3: Ne <<<")
lasso_results$Ne <- perform_lasso_detailed("Ne", c(abiotic_available,"aggregation_index"), complete_data)


#pglmm
sg_pglmm <- phyr::pglmm(as.formula(paste("sel.gain", "~", paste(paste(subset(lasso_results$sel.gain$selected, lasso_results$sel.gain$selected != "aggregation_index"), collapse = " + "),"(1|strain)",sep="+"))),
                        family = "gaussian",
                        data = complete_data, cov_ranef = list(strain = tree), bayes = FALSE)
summary(sg_pglmm)

sur_pglmm <- phyr::pglmm(as.formula(paste("SUR", "~", paste(paste(subset(lasso_results$SUR$selected, lasso_results$SUR$selected != "aggregation_index"), collapse = " + "),"(1|strain)",sep="+"))),
                         family = "gaussian",
                         data = complete_data, cov_ranef = list(strain = tree), bayes = FALSE)
summary(sur_pglmm)

ne_pglmm <- phyr::pglmm(as.formula(paste("Ne", "~", paste(paste(subset(lasso_results$Ne$selected, lasso_results$Ne$selected != "aggregation_index"), collapse = " + "),"(1|strain)",sep="+"))), 
                        family = "gaussian",
                        data = complete_data, cov_ranef = list(strain = tree), bayes = FALSE)
summary(ne_pglmm)


# Manual optimization
psem_model <- psem(lm(sel.gain ~ N + pH + SOC, complete_data), 
                   lm(SUR ~ elevation + precipitation + WHC + N + pH + sel.gain, complete_data), 
                   lm(Ne ~  precipitation + WHC + pH + SUR + sel.gain, complete_data)
)
summary(psem_model)


