\name{merit_de}
\alias{merit_de}
\title{MERIT-DE differential expression}
\description{Robust, FDR-safe, covariate-aware differential expression for bulk
RNA-seq: Winsorized robust transform + empirical-Bayes moderation + guarded
empirical-null calibration.}
\usage{merit_de(counts, group, covariates = NULL, transform = c("winsor","rint"),
  winsor_C = 2.5, mode = c("guard","base","full"), N_MIN = 20, TRIGGER = 1.5,
  SIGMA_CAP = 1.5, filtered = FALSE)}
\arguments{
\item{counts}{integer gene-by-sample count matrix (rows = genes).}
\item{group}{factor/vector of length ncol(counts) with exactly two levels.}
\item{covariates}{optional data.frame/vector of confounders to adjust for.}
\item{covariates}{optional data.frame/vector of confounders to adjust fr.}
#\item{transform}{"winsor" (default, robust + powerful) or "rint" (rank-INT).}
#\item{winsor_C}{Winsorization width in MAD units (default 2.5).}
#\item{mode}{"guard" (default), "base" (no empirical null), or "full".}
#\item{N_MIN}{minimum total sample size for empirical-null calibration (default 20).}
#\item{TRIGGER}{calibrate only if estimated sigma0 exceeds this (default 1.5).}
#\item{SIGMA_CAP}{cap on the empirical-null sigma0 (default 1.5).}
#\item{filtered}{set TRUE if `counts` is already gene-filtered.}
}
\value{data.frame with gene, stat, pvalue, padj.}
\examples{\dontrun{merit_de(counts, group)}}
