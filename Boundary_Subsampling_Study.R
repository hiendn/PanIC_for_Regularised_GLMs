#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "Boundary_Subsampling_Study.R"
replication_dir <- normalizePath(dirname(script_path))
source(file.path(replication_dir, "Auxiliary_Config.R"))
results_dir <- file.path(replication_dir, "results")
generated_dir <- results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

boundary_samples <- function(n, n_rep, seed, chunk_size = 2000L) {
  set.seed(seed)
  out <- numeric(n_rep)
  active <- logical(n_rep)
  cursor <- 1L
  while (cursor <= n_rep) {
    size <- min(chunk_size, n_rep - cursor + 1L)
    x <- matrix(rnorm(size * n), nrow = size, ncol = n)
    epsilon <- matrix(rnorm(size * n), nrow = size, ncol = n)
    y <- x + epsilon
    ## The proposition's scalar illustration has no intercept, so use the
    ## exact no-intercept least-squares estimator.
    beta_unconstrained <- rowSums(x * y) / rowSums(x^2)
    index <- cursor:(cursor + size - 1L)
    out[index] <- sqrt(n) * (pmin(beta_unconstrained, 1) - 1)
    active[index] <- beta_unconstrained >= 1
    cursor <- cursor + size
  }
  list(z = out, active = active)
}

fit_scalar_boundary <- function(x, y, bound = 1) {
  beta_unconstrained <- sum(x * y) / sum(x^2)
  min(beta_unconstrained, bound)
}

subsampling_samples <- function(x, y, b, n_sub, seed) {
  set.seed(seed)
  n <- length(x)
  beta_full <- fit_scalar_boundary(x, y)
  z <- numeric(n_sub)
  for (j in seq_len(n_sub)) {
    index <- sample.int(n, b, replace = FALSE)
    beta_sub <- fit_scalar_boundary(x[index], y[index])
    z[j] <- sqrt(b) * (beta_sub - beta_full)
  }
  list(z = z, beta_full = beta_full)
}

n_values <- c(100L, 500L)
empirical <- lapply(n_values, function(n) {
  boundary_samples(
    n, CONFIG$boundary_replications,
    seed = CONFIG$master_seed + 700000L + n
  )
})
boundary_summary <- do.call(rbind, lapply(seq_along(n_values), function(i) {
  z <- empirical[[i]]$z
  data.frame(
    n = n_values[i],
    replications = length(z),
    empirical_atom = mean(empirical[[i]]$active),
    atom_mcse = sqrt(mean(empirical[[i]]$active) *
                       (1 - mean(empirical[[i]]$active)) / length(z)),
    theoretical_atom = 0.5,
    q10 = unname(quantile(z, 0.10, type = 1)),
    q25 = unname(quantile(z, 0.25, type = 1)),
    theoretical_q10 = qnorm(0.10),
    theoretical_q25 = qnorm(0.25)
  )
}))
write.csv(boundary_summary, file.path(results_dir, "boundary_summary.csv"), row.names = FALSE)

n_full <- 2000L
b_values <- c(100L, 200L, 400L)
set.seed(CONFIG$master_seed + 800002L)
subsampling_x <- rnorm(n_full)
subsampling_y <- subsampling_x + rnorm(n_full)
subsampling <- lapply(seq_along(b_values), function(i) {
  subsampling_samples(
    subsampling_x, subsampling_y, b_values[i], CONFIG$subsampling_replications,
    seed = CONFIG$master_seed + 810000L + i
  )
})
cdf_points <- c(-1, -0.5, -0.1, 0.1)
subsampling_summary <- do.call(rbind, lapply(seq_along(b_values), function(i) {
  z <- subsampling[[i]]$z
  data.frame(
    n = n_full,
    b = b_values[i],
    b_over_n = b_values[i] / n_full,
    random_subsamples = length(z),
    beta_full = subsampling[[i]]$beta_full,
    atom_at_zero = mean(abs(z) < 1e-12),
    q10 = unname(quantile(z, 0.10, type = 1)),
    q25 = unname(quantile(z, 0.25, type = 1)),
    cdf_minus_1 = mean(z <= cdf_points[1]),
    cdf_minus_0_5 = mean(z <= cdf_points[2]),
    cdf_minus_0_1 = mean(z <= cdf_points[3]),
    cdf_plus_0_1 = mean(z <= cdf_points[4]),
    theoretical_cdf_minus_1 = pnorm(cdf_points[1]),
    theoretical_cdf_minus_0_5 = pnorm(cdf_points[2]),
    theoretical_cdf_minus_0_1 = pnorm(cdf_points[3]),
    theoretical_cdf_plus_0_1 = 1
  )
}))
write.csv(
  subsampling_summary,
  file.path(results_dir, "subsampling_summary.csv"), row.names = FALSE
)

figure_file <- file.path(generated_dir, "figure_boundary_subsampling.pdf")
pdf(figure_file, width = 10.2, height = 3.35, family = "Helvetica", pointsize = 9)
par(mfrow = c(1, 3), mar = c(3.6, 3.7, 2.0, 0.7), mgp = c(2.2, 0.7, 0), tcl = -0.25)
for (i in seq_along(n_values)) {
  z <- empirical[[i]]$z
  negative <- z[z < 0]
  breaks <- seq(-3.5, 0, length.out = 36L)
  negative_for_plot <- pmax(negative, breaks[1L] + 1e-10)
  h <- hist(negative_for_plot, breaks = breaks, plot = FALSE)
  plot(
    NA, xlim = c(-3.5, 0.12), ylim = c(0, 0.62),
    main = paste0("n = ", n_values[i]),
    xlab = expression(sqrt(n) * (hat(beta)[n] - beta^o)),
    ylab = if (i == 1L) "density / point mass" else ""
  )
  ## The bars integrate to the empirical negative-component probability rather
  ## than to one, leaving the remaining probability visible as the atom.
  rect(
    h$breaks[-length(h$breaks)], 0, h$breaks[-1L],
    h$density * length(negative) / length(z),
    col = "grey88", border = "grey60"
  )
  grid <- seq(-3.5, 0, length.out = 400L)
  lines(grid, dnorm(grid), lwd = 1.2)
  atom <- mean(empirical[[i]]$active)
  segments(0, 0, 0, atom, col = "#B2182B", lwd = 1.6)
  points(0, atom, pch = 16, cex = 0.7, col = "#B2182B")
  legend(
    "topleft", bty = "n", cex = 0.90,
    legend = c("Gaussian component", sprintf("atom = %.3f", atom)),
    lty = c(1, 1), pch = c(NA, 16),
    col = c("black", "#B2182B")
  )
}

plot(
  NA, xlim = c(-3.2, 0.35), ylim = c(0, 1.02),
  xlab = expression(x), ylab = "subsampling cdf", main = "n = 2000"
)
grid <- seq(-3.2, 0.35, length.out = 500L)
limit_cdf <- ifelse(grid < 0, pnorm(grid), 1)
lines(grid, limit_cdf, lwd = 1.6, col = "black")
cols <- c("#2166AC", "#67A9CF", "#D6604D")
for (i in seq_along(b_values)) {
  ez <- ecdf(subsampling[[i]]$z)
  lines(grid, ez(grid), col = cols[i], lwd = 1.1)
}
legend(
  "bottomright", bty = "n", cex = 0.90,
  legend = c("limit", paste0("b = ", b_values)),
  col = c("black", cols), lty = 1, lwd = c(1.6, rep(1.1, 3))
)
dev.off()
