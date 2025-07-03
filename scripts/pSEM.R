# Load required libraries
suppressPackageStartupMessages({
  library(piecewiseSEM)
  library(missForest)
  library(glmnet)
  library(broom)
  library(dplyr)
})

# Set your own working directory
#setwd("~")

cat("=== PIECEWISE SEM ANALYSIS - MICROBIAL EVOLUTION DATA ===\n")
cat("Enhanced version with detailed LASSO diagnostics and AI inclusion\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

# ===============================================================================
# DATA PREPARATION
# ===============================================================================

# Load and prepare data
data <- read.csv("SEM_data.csv", row.names = 1, stringsAsFactors = FALSE)
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
response_vars <- c("sel.gain", "SUR", "tajima.D")
abiotic_factors <- c("MAT", "precipitation", "elevation",
                     "evapo", "AI", "WHC", "N", "pH", "SOC")

# Prepare complete dataset
all_vars <- c(response_vars, abiotic_factors)
available_vars <- intersect(all_vars, colnames(data_standardized))
complete_data <- data_standardized[complete.cases(data_standardized[available_vars]), ]
abiotic_available <- intersect(abiotic_factors, colnames(complete_data))

cat("✓ Complete cases:", nrow(complete_data), "\n")
cat("✓ Available abiotic factors:", length(abiotic_available), "\n")

# Check AI availability
#if ("AI" %in% colnames(complete_data)) {
#  cat("✓ Aridity Index (AI) is available in the dataset\n")
#  ai_summary <- summary(complete_data$AI)
#  cat("AI summary statistics:\n")
#  print(ai_summary)
#} else {
#  cat("⚠ WARNING: AI not found in dataset!\n")
#}

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
  
  # Special AI analysis
  cat("\n--- ARIDITY INDEX (AI) ANALYSIS ---\n")
  if ("AI" %in% coef_detailed$Variable) {
    ai_row <- coef_detailed[coef_detailed$Variable == "AI", ]
    cat("AI coefficient at λ.min :", sprintf("%8.4f", ai_row$Coef_lambdamin), "\n")
    cat("AI coefficient at λ.1se :", sprintf("%8.4f", ai_row$Coef_lambda1se), "\n")
    
    # Format AI correlation with standard error
    if (!is.na(ai_row$Correlation_SE)) {
      cat("AI correlation with", response, ":", sprintf("%6.3f (±%5.3f)", ai_row$Correlation, ai_row$Correlation_SE), "\n")
    } else {
      cat("AI correlation with", response, ":", sprintf("%6.3f", ai_row$Correlation), "\n")
    }
    
    cat("AI selected at λ.min     :", ifelse(ai_row$Selected_lambdamin, "YES", "NO"), "\n")
    cat("AI selected at λ.1se     :", ifelse(ai_row$Selected_lambda1se, "YES ✓", "NO ✗"), "\n")
    
    if (!ai_row$Selected_lambda1se) {
      cat("⚠ AI NOT SELECTED by conservative λ.1se (will consider for post-hoc inclusion)\n")
    } else {
      cat("✓ AI SUCCESSFULLY SELECTED by λ.1se\n")
    }
  } else {
    cat("⚠ AI not in predictor set\n")
  }
  
  # Generate diagnostic plots
  png(paste0("lasso_outputs/detailed_", response, "_analysis.png"), 
      width = 16, height = 10, units = "in", res = 300)
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  
  # CV curve
  plot(cv_lasso, main = paste("CV Error:", response))
  abline(v = log(cv_lasso$lambda.1se), col = "red", lty = 2, lwd = 2)
  abline(v = log(cv_lasso$lambda.min), col = "blue", lty = 2, lwd = 2)
  legend("topright", c("λ.1se (conservative)", "λ.min"), col = c("red", "blue"), lty = 2, lwd = 2)
  
  # Coefficient paths
  plot(cv_lasso$glmnet.fit, xvar = "lambda", main = paste("Coefficient Paths:", response))
  abline(v = log(cv_lasso$lambda.1se), col = "red", lty = 2, lwd = 2)
  abline(v = log(cv_lasso$lambda.min), col = "blue", lty = 2, lwd = 2)
  
  # Selected variables
  if (length(selected_1se) > 0) {
    selected_coefs <- coef_detailed[coef_detailed$Selected_lambda1se & coef_detailed$Variable != "(Intercept)", ]
    barplot(abs(selected_coefs$Coef_lambda1se), names.arg = selected_coefs$Variable,
            main = "Selected Variables (λ.1se)", las = 2, cex.names = 0.8)
  } else {
    plot.new()
    text(0.5, 0.5, "No variables selected\nat λ.1se", cex = 2)
  }
  
  # AI coefficient tracking
  if ("AI" %in% coef_detailed$Variable && nrow(cv_lasso$glmnet.fit$beta) > 0) {
    if ("AI" %in% rownames(cv_lasso$glmnet.fit$beta)) {
      ai_path <- cv_lasso$glmnet.fit$beta["AI", ]
      plot(log(cv_lasso$glmnet.fit$lambda), ai_path, type = "l", lwd = 2,
           xlab = "log(λ)", ylab = "AI Coefficient", main = "AI Coefficient Path")
      abline(v = log(cv_lasso$lambda.1se), col = "red", lty = 2)
      abline(v = log(cv_lasso$lambda.min), col = "blue", lty = 2)
      abline(h = 0, col = "gray", lty = 3)
    }
  }
  
  # Correlation vs selection
  env_vars <- coef_detailed[coef_detailed$Variable != "(Intercept)", ]
  plot(abs(env_vars$Correlation), env_vars$Abs_coef_1se,
       xlab = "|Correlation|", ylab = "|Coefficient at λ.1se|",
       main = "Correlation vs Selection", pch = 19,
       col = ifelse(env_vars$Selected_lambda1se, "red", "gray"))
  if ("AI" %in% env_vars$Variable) {
    ai_point <- env_vars[env_vars$Variable == "AI", ]
    points(abs(ai_point$Correlation), ai_point$Abs_coef_1se, 
           pch = 19, cex = 2, col = "blue")
    text(abs(ai_point$Correlation), ai_point$Abs_coef_1se, 
         "AI", pos = 3, cex = 1.2, font = 2)
  }
  
  # Model complexity comparison
  barplot(c(length(selected_min), length(selected_1se)), 
          names.arg = c("λ.min", "λ.1se"),
          main = "Model Complexity", ylab = "Variables Selected")
  
  dev.off()
  
  # Export results
  write.csv(coef_detailed, paste0("lasso_outputs/", response, "_detailed_coefficients.csv"), row.names = FALSE)
  
  cat("\n--- SELECTION SUMMARY ---\n")
  cat("CONSERVATIVE SELECTION (λ.1se):", length(selected_1se), "variables\n")
  if (length(selected_1se) > 0) {
    cat("Selected:", paste(selected_1se, collapse = ", "), "\n")
  } else {
    cat("No variables selected (intercept-only)\n")
  }
  
  cat("Files generated: detailed_", response, "_analysis.png, ", response, "_detailed_coefficients.csv\n", sep = "")
  
  return(list(selected = selected_1se, cv_model = cv_lasso, coefficients = coef_detailed))
}

# ===============================================================================
# HIERARCHICAL LASSO SELECTION
# ===============================================================================

lasso_results <- list()

# Step 1: sel.gain
cat("\n>>> STEP 1: sel.gain <<<")
lasso_results$sel.gain <- perform_lasso_detailed("sel.gain", abiotic_available, complete_data)

# Step 2: SUR
cat("\n>>> STEP 2: SUR <<<")
lasso_results$SUR <- perform_lasso_detailed("SUR", abiotic_available, complete_data)

# Step 3: tajima.D
cat("\n>>> STEP 3: tajima.D <<<")
lasso_results$tajima.D <- perform_lasso_detailed("tajima.D", abiotic_available, complete_data)

# Extract selections
selected_vars <- list(
  sel.gain = lasso_results$sel.gain$selected,
  SUR = lasso_results$SUR$selected,
  tajima.D = lasso_results$tajima.D$selected
)

# ===============================================================================
# MODEL CONSTRUCTION
# ===============================================================================

cat("\n=== MODEL CONSTRUCTION ===\n")

build_model <- function(response, selected_predictors, fallback_predictors) {
  valid_predictors <- intersect(selected_predictors, colnames(complete_data))
  
  if (length(valid_predictors) == 0) {
    valid_predictors <- intersect(fallback_predictors, colnames(complete_data))
    cat("Using fallback for", response, "\n")
  }
  
  if (length(valid_predictors) > 0) {
    formula_str <- paste(response, "~", paste(valid_predictors, collapse = " + "))
    cat(formula_str, "\n")
    return(lm(as.formula(formula_str), data = complete_data))
  }
  return(NULL)
}

# Build initial models
models <- list(
  build_model("sel.gain", selected_vars$sel.gain),
  build_model("SUR", selected_vars$SUR),
  build_model("tajima.D", selected_vars$tajima.D)
)

# Remove NULL models
valid_models <- !sapply(models, is.null)
models <- models[valid_models]
model_names <- c("sel.gain", "SUR", "tajima.D")[valid_models]

cat("✓ Built", length(models), "models\n")

# ===============================================================================
# POST-LASSO AI INCLUSION (alternative)
# ===============================================================================

#cat("\n=== POST-LASSO AI INCLUSION ===\n")
#cat(paste(rep("=", 70), collapse = ""), "\n")

#force_ai_inclusion <- function(models, model_names, complete_data) {
#  ai_log <- data.frame(
#    Model = character(0),
#    Original_AI = logical(0),
#    AI_univar_coef = numeric(0),
#    AI_univar_p = numeric(0),
#    AI_added = logical(0),
#    Final_AI_coef = numeric(0),
#    Final_AI_p = numeric(0),
#    Justification = character(0),
#    stringsAsFactors = FALSE
#  )

#  modified_models <- models

#  for (i in 1:length(models)) {
#    if (!is.null(models[[i]])) {
#      model_name <- model_names[i]
#     current_terms <- attr(terms(models[[i]]), "term.labels")
#      ai_original <- "AI" %in% current_terms

#      cat("\n--- MODEL:", toupper(model_name), "---\n")
#      cat("Original predictors:", paste(current_terms, collapse = ", "), "\n")
#      cat("AI originally included:", ifelse(ai_original, "YES", "NO"), "\n")

# Test AI univariate
#      response_var <- all.vars(formula(models[[i]]))[1]
#      ai_univariate <- lm(formula(paste(response_var, "~ AI")), data = complete_data)
#      ai_p <- summary(ai_univariate)$coefficients["AI", "Pr(>|t|)"]
#      ai_coef <- summary(ai_univariate)$coefficients["AI", "Estimate"]

#      cat("AI univariate: coef =", sprintf("%.4f", ai_coef), ", p =", sprintf("%.4f", ai_p), "\n")

#      ai_added <- FALSE
#      final_coef <- NA
#      final_p <- NA
#      justification <- ""

#      if (!ai_original) {
#        if (ai_p < 0.2) {  # Ecological significance threshold
# Add AI
#          new_predictors <- c(current_terms, "AI")
#          new_formula <- as.formula(paste(response_var, "~", paste(new_predictors, collapse = " + ")))
#         modified_models[[i]] <- lm(new_formula, data = complete_data)
#        ai_added <- TRUE

#          final_summary <- summary(modified_models[[i]])
#          final_coef <- final_summary$coefficients["AI", "Estimate"]
#          final_p <- final_summary$coefficients["AI", "Pr(>|t|)"]
#          justification <- paste("Added: p =", sprintf("%.4f", ai_p), "< 0.2")

#          cat("✓ AI ADDED: final coef =", sprintf("%.4f", final_coef), ", p =", sprintf("%.4f", final_p), "\n")
#        } else {
#          justification <- paste("Not added: p =", sprintf("%.4f", ai_p), ">= 0.2")
#          cat("✗ AI NOT added:", justification, "\n")
#        }
#      } else {
#        final_summary <- summary(models[[i]])
#        final_coef <- final_summary$coefficients["AI", "Estimate"]
#        final_p <- final_summary$coefficients["AI", "Pr(>|t|)"]
#        justification <- "Originally selected by LASSO"
#        cat("✓ AI already included\n")
#      }

# Log results
#      ai_log <- rbind(ai_log, data.frame(
#        Model = model_name,
#        Original_AI = ai_original,
#        AI_univar_coef = ai_coef,
#        AI_univar_p = ai_p,
#        AI_added = ai_added,
#        Final_AI_coef = final_coef,
#        Final_AI_p = final_p,
#        Justification = justification,
#        stringsAsFactors = FALSE
#      ))
#   }
#}

# Summary
#cat("\n", paste(rep("-", 70), collapse = ""), "\n")
#cat("AI INCLUSION SUMMARY\n")
#print(ai_log)

#  ai_final <- sum(!is.na(ai_log$Final_AI_coef))
#  ai_sig <- sum(ai_log$Final_AI_p < 0.05, na.rm = TRUE)

#  cat("\nFINAL STATUS:\n")
#  cat("- Models with AI included:", ai_final, "\n")
# cat("- Models with significant AI:", ai_sig, "\n")

#  if (ai_final > 0) {
#   cat("✓ SUCCESS: AI included in", ai_final, "model(s)\n")
#} else {
# cat("⚠ WARNING: AI not included in any model\n")
#}

#write.csv(ai_log, "AI_inclusion_log.csv", row.names = FALSE)
#cat("✓ Log exported: AI_inclusion_log.csv\n")

#return(list(models = modified_models, log = ai_log))
#}

# Apply AI inclusion
#ai_results <- force_ai_inclusion(models, model_names, complete_data)
#models <- ai_results$models

# ===============================================================================
# PIECEWISE SEM ANALYSIS
# ===============================================================================

cat("\n=== PIECEWISE SEM ANALYSIS ===\n")

build_psem <- function(models) {
  if (length(models) == 3) {
    return(psem(models[[1]], models[[2]], models[[3]]))
  } else if (length(models) == 2) {
    return(psem(models[[1]], models[[2]]))
  } else {
    return(psem(models[[1]]))
  }
}

tryCatch({
  psem_model <- build_psem(models)
  psem_summary <- summary(psem_model)
  
  cat("✓ pSEM created\n")
  cat("Fisher's C =", round(psem_summary$Cstat$Fisher.C, 3), 
      "(P =", round(psem_summary$Cstat$P.Value, 4), ")\n")
  cat("AIC =", round(psem_summary$AIC$AIC, 2), "\n")
  
  # Check AI paths
  #if ("coefficients" %in% names(psem_summary)) {
  #  ai_paths <- psem_summary$coefficients[grepl("AI", psem_summary$coefficients$Predictor), ]
  #  if (nrow(ai_paths) > 0) {
  #    cat("\n✓ AI PATH COEFFICIENTS:\n")
  #    print(ai_paths[, c("Response", "Predictor", "Estimate", "Std.Error", "P.Value")])
  #  }
  #}
  
  # Export results
  capture.output({
    cat("===============================================================================\n")
    cat("PIECEWISE SEM ANALYSIS - WITH AI INCLUSION\n")
    cat("===============================================================================\n\n")
    cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
    cat("Sample Size:", nrow(complete_data), "\n\n")
    
    cat("\nFINAL MODELS:\n")
    for (i in 1:length(models)) {
      cat("Model", i, ":", paste(deparse(formula(models[[i]])), collapse = ""), "\n")
    }
    
    cat("\nPSEM RESULTS:\n")
    print(psem_summary)
    
    cat("\nINDIVIDUAL MODEL SUMMARIES:\n")
    for (i in 1:length(models)) {
      cat("\n--- MODEL", i, "---\n")
      print(summary(models[[i]]))
    }
  }, file = "pSEM_RESULTS_with_AI.txt")
  
  save.image("pSEM_WORKSPACE_with_AI.RData")
  
  cat("\n✓ Results exported:\n")
  cat("  - pSEM_RESULTS_with_AI.txt\n")
  cat("  - lasso_outputs/ (detailed diagnostics)\n")
  cat("  - pSEM_WORKSPACE_with_AI.RData\n")
  
}, error = function(e) {
  cat("Error:", e$message, "\n")
  save.image("pSEM_ERROR.RData")
})

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("COMPLETE pSEM ANALYSIS WITH AI INCLUSION - FINISHED\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("✓ Conservative λ.1se selection maintained\n")
cat("✓ Detailed LASSO diagnostics generated\n")
cat("✓ pSEM analysis completed\n") 



# Manual optimization
summary(psem_model)
psem_model <- psem(lm(sel.gain ~ N + pH + SOC, complete_data), 
                   lm(SUR ~ elevation + evapo + AI + WHC + N + pH + sel.gain, complete_data), 
                   lm(tajima.D ~ AI + N + pH + SUR + sel.gain, complete_data)
)


