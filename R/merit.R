#' Estimate the empirical null (Efron 2004 central matching)
#'
#' Fits a Gaussian to the central part of the genome-wide z statistics to recover
#' the true null N(delta0, sigma0^2) when the theoretical N(0,1) is miscalibrated.
#' @param z numeric vector of z statistics
#' @param p_central central fraction used for matching (default 0.5)
#' @return list(delta0, sigma0)
#' @export
emp_null_cm <- function(z, p_central = 0.5){
  z <- z[is.finite(z)]
  qs <- stats::quantile(z, c((1-p_central)/2, 1-(1-p_central)/2), names = FALSE)
  res <- tryCatch({
    d <- stats::density(z, n = 1024)
    sel <- d$x >= qs[1] & d$x <= qs[2] & d$y > 0
    if (sum(sel) >= 10){
      x <- d$x[sel]; ly <- log(d$y[sel])
      fit <- stats::lm(ly ~ x + I(x^2)); b1 <- coef(fit)[2]; b2 <- coef(fit)[3]
      if (is.finite(b2) && b2 < 0){
        sigma0 <- sqrt(-1/(2*b2)); list(delta0 = as.numeric(b1*sigma0^2), sigma0 = as.numeric(sigma0))
      } else NULL
    } else NULL
  }, error = function(e) NULL)
  if (is.null(res)){
    zc <- z[z>=qs[1] & z<=qs[2]]
    res <- list(delta0 = stats::median(z), sigma0 = max(stats::mad(zc), 1e-6))
  }
  res$sigma0 <- max(res$sigma0, 0.5); res
}

#' MERIT-DE differential expression
#'
#' @param counts integer gene-by-sample count matrix (rows = genes).
#' @param group factor/vector of length ncol(counts) with exactly two levels.
#' @param covariates optional data.frame/vector of confounders to adjust for.
#' @param transform "winsor" (default, robust + powerful) or "rint" (rank-INT).
#' @param winsor_C Winsorization width in MAD units (default 2.5).
#' @param mode "guard" (default), "base" (no empirical null), or "full".
#' @param N_MIN minimum total sample size for empirical-null calibration (20).
#' @param TRIGGER calibrate only if estimated sigma0 exceeds this (default 1.5).
#' @param SIGMA_CAP cap on the empirical-null sigma0 (default 1.5).
#' @param min_count filterByExpr is applied unless filtered = TRUE.
#' @param filtered set TRUE if `counts` is already gene-filtered.
#' @return data.frame with gene, stat (z), pvalue, padj (BH).
#' @examples
#' \dontrun{ res <- merit_de(counts, group); head(res[order(res$padj),]) }
#' @export
merit_de <- function(counts, group, covariates = NULL,
                     transform = c("winsor","rint"), winsor_C = 2.5,
                     mode = c("guard","base","full"), N_MIN = 20,
                     TRIGGER = 1.5, SIGMA_CAP = 1.5, filtered = FALSE){
  transform <- match.arg(transform); mode <- match.arg(mode)
  group <- factor(group)
  if (nlevels(group) != 2) stop("group must have exactly two levels")
  y <- edgeR::DGEList(counts, group = group)
  if (!filtered) y <- y[edgeR::filterByExpr(y, group = group), , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  genes <- rownames(logcpm); n <- ncol(logcpm)
  ## robust transform
  if (transform == "rint"){
    X <- t(apply(logcpm, 1, function(x){ qnorm((rank(x, ties.method="average")-0.5)/n) }))
  } else {
    X <- t(apply(logcpm, 1, function(x){
      md <- stats::median(x); s <- stats::mad(x); if (s < 1e-6) s <- stats::sd(x)+1e-6
      pmin(pmax(x, md-winsor_C*s), md+winsor_C*s) }))
  }
  design <- if (is.null(covariates)) stats::model.matrix(~ group) else
            stats::model.matrix(~ group + ., data = data.frame(covariates))
  fit <- limma::eBayes(limma::lmFit(X, design))
  tstat <- fit$t[,2]; df <- fit$df.total
  z <- stats::qnorm(stats::pt(tstat, df = df))
  z[!is.finite(z)] <- sign(tstat[!is.finite(z)]) * 8
  en <- emp_null_cm(z)
  if (mode == "base"){ d0 <- 0; s0 <- 1 }
  else if (mode == "full"){ d0 <- en$delta0; s0 <- en$sigma0 }
  else { reliable <- n >= N_MIN
    s0 <- if (reliable && en$sigma0 > TRIGGER) min(en$sigma0, SIGMA_CAP) else 1
    d0 <- if (reliable && abs(en$delta0) > 0.1) max(min(en$delta0, 0.5), -0.5) else 0 }
  zc <- (z - d0)/s0
  p <- 2*stats::pnorm(-abs(zc))
  data.frame(gene = genes, stat = zc, pvalue = p,
             padj = stats::p.adjust(p, "BH"), row.names = NULL)
}
