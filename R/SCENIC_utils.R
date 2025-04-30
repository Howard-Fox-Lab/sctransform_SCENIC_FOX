#' Binomial variance-stabilizing transformation for AUCell scores
#' @generated text
#' @param auc_scores A numeric vector or matrix of AUCell scores in [0, 1]
#' @return A transformed vector or matrix using the arcsin(sqrt(p)) transform
#'
#' @section Details:
#' AUCell scores are bounded between 0 and 1 and reflect enrichment probabilities
#' of regulon activation across cells. These scores exhibit heteroskedastic variance,
#' especially near the boundaries of 0 and 1. This function applies the binomial
#' variance-stabilizing transformation (VST), which stabilizes variance across the
#' dynamic range of AUCell outputs.
#'
#' This approach is analogous to the SCTransform method for gene counts but adapted
#' for regulon activation scores, which behave like soft Bernoulli trials.
#'
#' This transformation is essential for PCA, batch correction, and downstream
#' distance-based analyses on AUCell-derived matrices.
#'
#' @importFrom magrittr %>%
#' @importFrom stats asin sqrt

library(AUCell)
library(ggplot2)

## Helper Functions

remove_extended <- function(df) {
  df = df[!grepl("extended", rownames(df)), ]
}

bernoolize <- function(x) {
  stopifnot(all(x >= 0 & x <= 1))
  asin(sqrt(x))
}

bern_var <- function(x) {
  p <- colMeans(x)
  p * (1 - p)
}

pearson_resid_auc <- function(p_obs, p_pred) {
  stopifnot(length(p_obs) == length(p_pred))
  (p_obs - p_pred) / sqrt(p_obs * (1 - p_obs) + 1e-6)
}

fit_regulon_model <- function(df, formula = value ~ regulon_size) {
  lm(formula, data = df)
}

vlnplot <- function(df, regulon_name) {
  ggplot(df, aes(x = label, y = .data[[regulon_name]])) +
    geom_violin(trim = FALSE, fill = "lightblue") +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_minimal() +
    labs(title = gsub(" \\(.*\\)", "", regulon_name),  # Clean up title, remove (3214g)
         x = "Cell Type", y = "AUC Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

