# Install and load plotly
install.packages("plotly")
library(plotly)

# Example matrix
set.seed(42)
mat <- matrix(rnorm(100), nrow = 10)

# Create a 3D bar plot
plot_ly(z = mat, type = "surface", colorscale = "Viridis") %>%
  layout(title = "3D Barplot of a Matrix")

gen1 <- spouse_pairs(pop, n_iter = 1e5)

gen1_df <- data.frame(names(gen1) |> as.integer(), gen1,
                      pop$get_sex(0)$phi1, pop$get_sex(1)$phi1)

colnames(gen1_df) <- c("m", "f", "phi_m", "phi_f")

plot(gen1_df$phi_m, gen1_df$phi_f)
ggplot(gen1_df, aes(x = phi_m, y = phi_f)) + geom_point(position = "jitter") +
  geom_smooth(method = "lm")
table(gen1_df$phi_m, gen1_df$phi_f) |> prop.table()
