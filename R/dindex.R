# dindex.R -- two-group log-ratio index D = (P - M) / (P + M) = tanh(0.5 * log(P / M))
# MIT License (see LICENSE). Base R only; the phyloseq helper is used if phyloseq is installed.
# Source with: source("dindex.R")

d_index <- function(plus, minus) {
  plus <- as.numeric(plus); minus <- as.numeric(minus)
  if (any(plus < 0, na.rm = TRUE) || any(minus < 0, na.rm = TRUE)) stop("abundances must be non-negative")
  den <- plus + minus
  ifelse(den > 0, (plus - minus) / ifelse(den > 0, den, 1), NA_real_)
}

log_ratio <- function(plus, minus, pseudocount = NULL) {
  plus <- as.numeric(plus); minus <- as.numeric(minus)
  if (any(plus < 0, na.rm = TRUE) || any(minus < 0, na.rm = TRUE)) stop("abundances must be non-negative")
  if (!is.null(pseudocount) && (!is.finite(pseudocount) || pseudocount <= 0)) stop("pseudocount must be a finite positive number")
  if (!is.null(pseudocount)) { plus <- plus + pseudocount; minus <- minus + pseudocount }
  log(plus) - log(minus)
}

d_from_log_ratio <- function(lr) tanh(lr / 2)

# abundance: samples x features matrix/data.frame (features_as = "columns") or features x samples ("rows")
# taxonomy : named character vector, names = feature IDs, values = group labels (e.g. order)
d_from_table <- function(abundance, taxonomy, plus, minus, features_as = c("columns", "rows"),
                         min_share = NULL, min_total = NULL, pseudocount = NULL) {
  features_as <- match.arg(features_as)
  if (!is.null(min_share) && (min_share < 0 || min_share > 1)) stop("min_share must be between 0 and 1")
  if (!is.null(min_total) && min_total < 0) stop("min_total must be non-negative")
  tab <- as.matrix(abundance)
  if (features_as == "rows") tab <- t(tab)
  if (length(intersect(plus, minus)) > 0) stop("plus and minus groups overlap")
  mapped <- intersect(colnames(tab), names(taxonomy))
  if (length(mapped) == 0) stop("no abundance features found in taxonomy")
  grp <- taxonomy[mapped]
  missing <- setdiff(c(plus, minus), unique(grp))
  if (length(missing) > 0) stop("group labels not present in taxonomy: ", paste(missing, collapse = ", "))
  vals <- tab[, mapped, drop = FALSE]
  if (any(vals < 0, na.rm = TRUE)) stop("abundances must be non-negative")
  P <- rowSums(vals[, grp %in% plus, drop = FALSE]); M <- rowSums(vals[, grp %in% minus, drop = FALSE])
  total <- rowSums(vals)
  out <- data.frame(sample = rownames(tab), plus = P, minus = M, total = total, D = d_index(P, M),
                    log_ratio = log_ratio(P, M, pseudocount), share = ifelse(total > 0, (P + M) / total, NA_real_),
                    stringsAsFactors = FALSE)
  qc <- rep("ok", nrow(out))
  if (!is.null(min_total)) qc[(P + M) < min_total] <- "low_total"
  if (!is.null(min_share)) qc[qc == "ok" & (is.na(out$share) | out$share < min_share)] <- "low_share"
  qc[(P + M) == 0] <- "undefined"
  out$qc <- qc
  rownames(out) <- out$sample
  out
}

# Beta (Jeffreys) credible interval for D from read counts: p = P/(P+M) ~ Beta(nP + prior, nM + prior), D = 2p - 1
d_ci_counts <- function(n_plus, n_minus, level = 0.95, prior = 0.5) {
  if (any(n_plus < 0, na.rm = TRUE) || any(n_minus < 0, na.rm = TRUE)) stop("counts must be non-negative")
  if (!is.finite(level) || level <= 0 || level >= 1) stop("level must be between 0 and 1")
  if (!is.finite(prior) || prior <= 0) stop("prior must be a finite positive number")
  a <- n_plus + prior; b <- n_minus + prior; alpha <- (1 - level) / 2
  data.frame(D = d_index(n_plus, n_minus), D_lower = 2 * qbeta(alpha, a, b) - 1, D_upper = 2 * qbeta(1 - alpha, a, b) - 1)
}

# phyloseq helper: sums taxa by a tax_table rank
d_from_phyloseq <- function(physeq, rank, plus, minus, ...) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("phyloseq is not installed")
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (phyloseq::taxa_are_rows(physeq)) otu <- t(otu)
  tax <- as.character(phyloseq::tax_table(physeq)[, rank]); names(tax) <- phyloseq::taxa_names(physeq)
  d_from_table(otu, tax[!is.na(tax)], plus, minus, ...)
}

# Driver value at which a linear fit of D crosses zero, with a bootstrap 95% interval (optionally resampling groups)
d_zero_point <- function(data, driver, d_col = "D", n_boot = 2000, group = NULL, seed = 0) {
  cols <- c(d_col, driver, if (!is.null(group)) group)
  df <- data[stats::complete.cases(data[, cols, drop = FALSE]), , drop = FALSE]
  if (n_boot < 2 || n_boot != as.integer(n_boot)) stop("n_boot must be an integer of at least 2")
  if (nrow(df) < 3 || length(unique(df[[driver]])) < 2) stop("at least three complete observations and two driver values are required")
  if (!is.null(group) && length(unique(df[[group]])) < 2) stop("group bootstrap requires at least two groups")
  x0 <- function(d) { b <- stats::coef(stats::lm(d[[d_col]] ~ d[[driver]])); -b[1] / b[2] }
  set.seed(seed)
  boots <- replicate(n_boot, {
    if (is.null(group)) x0(df[sample.int(nrow(df), replace = TRUE), ])
    else { u <- unique(df[[group]]); pick <- sample(u, length(u), replace = TRUE)
           x0(do.call(rbind, lapply(pick, function(g) df[df[[group]] == g, ]))) }
  })
  c(zero_point = unname(x0(df)), ci_low = unname(stats::quantile(boots, 0.025, na.rm = TRUE)),
    ci_high = unname(stats::quantile(boots, 0.975, na.rm = TRUE)), n = nrow(df))
}
