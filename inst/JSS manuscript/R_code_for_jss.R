# Preliminaries
library(CPN)  # Assuming your package is installed

# Simulate data
set.seed(123)
data <- simulate_cpn_data()
head(data)

# Basic visualization
stripchart(
  y ~ x1, data = data,
  vertical = TRUE, method = "jitter",
  pch = 19, cex = 0.6,
  col = c("blue", "red"),
  xlab = "x1", ylab = "y",
  main = "Scatter Plot of y by Group x1"
)

# Fit CPN model
fit <- cpn(y ~ x1 + x2, data = data)

# Model summary
summary(fit)

# Extract coefficients
coef(fit)               # Full coefficients
coef(fit, full = FALSE) # Regression coefficients only

# Diagnostics plots
plot(fit)

# Variance-covariance matrix
vcov(fit)

# ANOVA (Type I deviance)
anova(fit)

# Predictions on new data
new_df <- data.frame(
  x1 = c("A", "B"),
  x2 = c(0.5, -0.3)
)
predict(fit, newdata = new_df, type = "response", interval = "confidence")

# Compare to classical linear model
lm_fit <- lm(y ~ x1 + x2, data = data)
anova_lm <- anova(lm_fit)
anova_cpn <- anova(fit)

# Extract p-values for x1 safely
print("CPN ANOVA row names:")
print(anova_cpn$Term)
print("LM ANOVA row names:")
print(rownames(anova_lm))

cpn_row <- grep("x1", anova_cpn$Term)
if (length(cpn_row) == 1) {
  cpn_pval <- anova_cpn[cpn_row, "Pr(>Chi)"]
} else {
  cpn_pval <- NA
}

lm_row <- grep("^x1", rownames(anova_lm))
if (length(lm_row) == 1) {
  lm_pval <- anova_lm[lm_row, "Pr(>F)"]
} else {
  lm_pval <- NA
}

# Format p-values
format_pval <- function(p) {
  if (is.na(p)) NA
  else if (p < 0.001) "<0.001"
  else sprintf("%.3f", p)
}

lm_pval_tex <- format_pval(lm_pval)
cpn_pval_tex <- format_pval(cpn_pval)

# Print p-values
cat("Linear model p-value for x1:", lm_pval_tex, "\n")
cat("CPN model p-value for x1:", cpn_pval_tex, "\n")
