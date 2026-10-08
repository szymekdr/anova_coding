## anova_coding_demo.R
## Companion to anova_coding_manual.html (Szymek Drobniak).
## Simulated nestling body mass in 4 habitats; shows how each coding scheme
## changes the coefficients but not the fit, and how B^{-1} gives their meaning.
## Note: R's RNG differs from the JavaScript generator in the HTML page,
## so the numbers will not match the page exactly; the algebra does.

set.seed(20261007)
habitats <- c("Deciduous", "Mixed", "Coniferous", "Urban park")
true_mu  <- c(15.4, 16.9, 14.2, 13.6)
n        <- c(10, 10, 10, 10)        # try c(4, 12, 18, 7) for an unbalanced design
sd_y     <- 1.15

d <- data.frame(
  habitat = factor(rep(habitats, n), levels = habitats),
  mass    = rnorm(sum(n), rep(true_mu, n), sd_y)
)

## ---- helper: fit with a given coding and report everything -----------------
show_coding <- function(dat, C, label) {
  contrasts(dat$habitat) <- C
  fit <- lm(mass ~ habitat, data = dat)
  B   <- cbind(Intercept = 1, C)            # one row per group: mu = B %*% beta
  Bi  <- solve(B)                           # hypothesis matrix: beta = Bi %*% mu
  ybar <- tapply(dat$mass, dat$habitat, mean)
  cat("\n==========", label, "==========\n")
  cat("Coding matrix C:\n");            print(MASS::fractions(C))
  cat("Hypothesis matrix B^{-1} (rows = weights on group means):\n")
  print(MASS::fractions(Bi))
  cat("lm() coefficients vs B^{-1} %*% group means:\n")
  print(cbind(lm = coef(fit), B_inv_ybar = as.vector(Bi %*% ybar)))
  cat("Fitted group means (invariant):", round(unique(fitted(fit)), 3), "\n")
  invisible(fit)
}

k <- nlevels(d$habitat)

## ---- reduced (full-rank) codings ------------------------------------------
fits <- list(
  treatment = show_coding(d, contr.treatment(k), "treatment (dummy, ref = level 1)"),
  sum       = show_coding(d, contr.sum(k),       "sum-to-zero (+1/0/-1)"),
  helmert   = show_coding(d, contr.helmert(k),   "Helmert"),
  sdif      = show_coding(d, MASS::contr.sdif(k),"successive differences"),
  poly      = show_coding(d, contr.poly(k),      "orthogonal polynomial")
)

## all fits identical: same residual SS and same omnibus F
sapply(fits, function(f) c(RSS = sum(resid(f)^2), F = anova(f)$`F value`[1]))

## ---- cell-means model ------------------------------------------------------
coef(lm(mass ~ 0 + habitat, data = d))      # = group means

## ---- effects model is rank-deficient ---------------------------------------
X_full <- cbind(1, model.matrix(~ 0 + habitat, data = d))
qr(X_full)$rank                             # k, not k + 1

## ---- two groups: the factor-of-2 trap ---------------------------------------
d2 <- droplevels(subset(d, habitat %in% habitats[1:2]))
contrasts(d2$habitat) <- contr.sum(2)       # Deciduous = +1, Mixed = -1
b_pm1  <- coef(lm(mass ~ habitat, d2))[2]
contrasts(d2$habitat) <- -contr.sum(2) / 2  # Deciduous = -1/2, Mixed = +1/2
b_half <- coef(lm(mass ~ habitat, d2))[2]
diff_means <- diff(tapply(d2$mass, d2$habitat, mean))
c(plus_minus_1 = b_pm1, plus_minus_half = b_half, mean_difference = diff_means)

## ---- unbalanced: intercept of centred codings != overall mean --------------
d_unb <- d[unlist(mapply(function(g, m) which(d$habitat == g)[seq_len(m)],
                         habitats, c(4, 10, 10, 7))), ]
contrasts(d_unb$habitat) <- contr.sum(k)
c(intercept_sum_coding = unname(coef(lm(mass ~ habitat, d_unb))[1]),
  mean_of_group_means  = mean(tapply(d_unb$mass, d_unb$habitat, mean)),
  overall_sample_mean  = mean(d_unb$mass))

## ---- sampling correlation of slope estimates --------------------------------
slope_cor <- function(C, n) {
  Bi <- solve(cbind(1, C))
  V  <- Bi %*% diag(1 / n) %*% t(Bi)
  cov2cor(V)[-1, -1]
}
round(slope_cor(contr.sum(3), c(10, 10, 10)), 3)      # correlated (not orthogonal)
round(slope_cor(contr.helmert(3), c(10, 10, 10)), 3)  # uncorrelated
