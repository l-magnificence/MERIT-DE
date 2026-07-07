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
\item{counts}{integer gene-by-sample count matrix.}
\item{group}{two-level factor of conditions.}
\item{covariates}{optional confounders to adjust for.}
}
\value{data.frame with gene, stat, pvalue, padj.}
\examples{\dontrun{merit_de(counts, group)}}
