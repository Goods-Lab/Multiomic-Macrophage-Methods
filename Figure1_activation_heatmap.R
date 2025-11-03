library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggplot2)
library(viridis)

#### Activation Heatmap with cells below LOD greyed out #########

# Read data
fileint <- '~Figure1_activation/activation.csv' #Include the appropriate file path
data_int <- read.csv(fileint, header=TRUE, row.names = 1)

# Read metadata
metadata <- '~Figure1_activation/activation_metadata.csv' #Include the appropriate file path
metadata <- read.csv(metadata)

# Extract relevant columns for annotation
annotation_cols <- metadata[, c("Sex", "Activation", "CellType")]

# Convert annotation_cols to data frame 
annotation_cols <- as.data.frame(annotation_cols)
annotation_cols$Sex <- as.factor(annotation_cols$Sex)
annotation_cols$Activation <- as.factor(annotation_cols$Activation)
annotation_cols$CellType <- as.factor(annotation_cols$CellType)

# Define LOD vector corresponding to each cytokine
LOD_vector <- c(91.03761706, 62.456, 93.258, 107.037, 149.607, 175.187, 206.116, 
                73.0536, 718.839)  # LOD for each cytokine

# Replace values below LOD with NA
data_int_na <- data_int
for (i in 1:nrow(data_int)) {
  data_int_na[i, data_int[i, ] < LOD_vector[i]] <- NA
}

# Scale the data by row, excluding NA values
data_int_matrix_scaled <- t(apply(data_int_na, 1, function(x) {
  non_na_values <- x[!is.na(x)]
  if (length(non_na_values) > 1) {
    scaled_values <- scale(non_na_values)
    x[!is.na(x)] <- scaled_values
  }
  return(x)
}))

# Define color palette from red to blue for z-score
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(25)

# Extract relevant columns for annotation
annotation_cols <- metadata[, c("Sex", "Activation")]

# Convert annotation_cols to factors
annotation_cols$Sex <- as.factor(annotation_cols$Sex)
annotation_cols$Activation <- as.factor(annotation_cols$Activation)

# Define color gradients for Sex and Activation
sex_colors <- colorRampPalette(c("pink2", "mediumseagreen"))(length(levels(annotation_cols$Sex)))
activation_colors <- colorRampPalette(c("grey","darkred", "skyblue3" ))(length(levels(annotation_cols$Activation)))

# Create named vectors for annotation colors
sex_col_colors <- setNames(sex_colors, levels(annotation_cols$Sex))
activation_col_colors <- setNames(activation_colors, levels(annotation_cols$Activation))

# Create heatmap annotation with borders
ha <- HeatmapAnnotation(
  df = annotation_cols,
  col = list(Sex = sex_col_colors, Activation = activation_col_colors),
  annotation_name_side = "left",
  gp = gpar(col = "black", lwd = 1) # Customize border color and line width
)

# Define a custom cell function to add borders
cell_border_function <- function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width = width, height = height,
            gp = gpar(col = "black", fill = NA, lwd = 0.5)) 
}

# Draw the heatmap with ComplexHeatmap and add annotations and legend for LOD
draw(Heatmap(data_int_matrix_scaled,
             name = "z-score",
             top_annotation = ha,
             cluster_columns = TRUE,
             cluster_rows= TRUE,
             #clustering_distance_rows = euclidean,
             clustering_method_rows = "complete",
             col = heatmap_colors,  # Use the defined red to blue color palette
             na_col = "grey",
             row_names_gp = gpar(fontsize = 10), 
             column_names_gp = gpar(fontsize = 10),
             cell_fun = cell_border_function, # Apply custom cell function
             heatmap_legend_param = list(
               title = "z-score",
               at = seq(-2, 2, by = 1)  # Adjust the legend breaks as needed
             )))