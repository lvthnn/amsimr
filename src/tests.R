# Install and load plotly
install.packages("plotly")
library(plotly)

# Example matrix
set.seed(42)
mat <- matrix(rnorm(100), nrow = 10)

# Create a 3D bar plot
plot_ly(z = mat, type = "surface", colorscale = "Viridis") %>%
  layout(title = "3D Barplot of a Matrix")