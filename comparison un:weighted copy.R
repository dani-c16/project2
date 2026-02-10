calc_mC_cov_unweighted(bs_obj, wholechromo1)
calc_mC_cov_weighted(bs_obj, wholechromo2)
my_promoters
calc_mC_cov_unweighted(bs_obj, my_promoters)
calc_mC_cov_weighted(bs_obj, my_promoters)
calc_mC_cov_perbase(bs_obj, my_promoters)
library(ggplot2)
ahhh <- wholechromo1[1:22]
wholechromo1M <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "M", what = "perRegionTotal")
wholechromo2M <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "M", what = "perRegionAverage")
wholechromo2C <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "Cov", what = "perRegionAverage")

library(ggplot2)

# Combine both results into a single dataframe for easier plottihist
hist(wholechromo1M[1:22], col ="purple")
plot(wholechromo2M[1:22], wholechromo2C[1:22], 
     xlab = "M", ylab = "Coverage", main = "Unweighted Methylation of Chromosomes 1 to 22", col = "purple")
model <- lm(wholechromo2C[1:22] ~ wholechromo2M[1:22])

# Add line of best fit
abline(model, col = "lightblue", lwd = 2)  # lwd controls the line width


# Sample data
set.seed(42)  # For reproducibility
x <- seq(-3, 3, length.out = 50) 
y <- x^3 - 2*x^2 + rnorm(50, 0, 2)  # Cubic relationship with noise

# Basic scatter plot
plot(x, y, xlab = "X Axis Title", ylab = "Y Axis Title", main = "Scatter Plot with Curve of Best Fit", col = "blue", pch = 19)

#obj for weighted data
coverage_df <- data.frame(
  Position = seq_along(wholechromo2M[1:22]),
  Coverage = as.numeric(wholechromo2M[1:22])
)
coverage_df2 <- data.frame(
  Position = seq_along(wholechromo2C[1:22]),
  Coverage = as.numeric(wholechromo2C[1:22])
)
# Load ggplot2 for visualization
library(ggplot2)

# Visualize the coverage
ggplot(coverage_df, aes(x = Position, y = Coverage)) +
  geom_line(color = "pink") +
  labs(title = "Unweighted Methylation Coverage Across Genome", 
       x = "Genomic Position", y = "Methylation Coverage") + theme_dark()
#graph comparing M value to Cov value for weighted
ggplot() +
  geom_line(data = coverage_df, aes(x = Position, y = Coverage), color = "pink", size = 1) + 
  geom_line(data = coverage_df2, aes(x = Position, y = Coverage), color = "purple", size = 1) +
  labs(title = "Methylation and Coverage of Chromosomes 1 to 22 (Weighted)", colour = "legend") +
  scale_color_manual(name = "Legend", values = c("Methylation" = "pink", "Coverage" = "purple")) +
  theme_dark()


legend("topleft", legend = c("Methylation", "Coverage"), col = c("pink", "purple"), lty = 1, bty = "n")





