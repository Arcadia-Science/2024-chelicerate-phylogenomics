require(grDevices)
require(Matrix)

# Some color palettes:
# palettes as R vectors
accent_v2 <-
  c(
    "#5088C5",
    "#F28360",
    "#3B9886",
    "#F7B846",
    "#7A77AB",
    "#F898AE",
    "#C6E7F4",
    "#F8C5C1",
    "#B5BEA4",
    "#F5E4BE",
    "#DCBFFC",
    "#F5CBE4"
  )

accent_v2_expanded <-
  c(
    "#97CD78",
    "#73B5E3",
    "#FFB984",
    "#BAB0A8",
    "#C85152",
    "#8A99AD",
    "#D1EADF",
    "#BABEE0",
    "#F1E8DA",
    "#DAD3C7",
    "#DA9085",
    "#B6C8D4"
  )

accent_ordered <-
  c(
    "#5088C5",
    "#F28360",
    "#F7B846",
    "#97CD78",
    "#7A77AB",
    "#F898AE",
    "#3B9886",
    "#C85152",
    "#73B5E3",
    "#BAB0A8",
    "#8A99AD",
    "#FFB984"
  )

light_ordered <-
  c(
    "#C6E7F4",
    "#F8C5C1",
    "#F5E4BE",
    "#B5BEA4",
    "#DCBFFC",
    "#B6C8D4",
    "#DAD3C7",
    "#DA9085",
    "#F5CBE4",
    "#D1EADF",
    "#BABEE0",
    "#F1E8DA"
  )

accent_all <-
  c(accent_ordered, light_ordered)
accent_v3 <-
  c(
    "#5088C5",
    "#C85152",
    "#F7B846",
    "#3B9886",
    "#7A77AB",
    "#73B5E3",
    "#F28360",
    "#FFB984",
    "#97CD78",
    "#DCBFFC",
    "#C6E7F4",
    "#F898AE",
    "#F5E4BE",
    "#D1EADF",
    "#BABEE0"
  )
accent_v4 <-
  c(
    "#5088C5",
    "#C85152",
    "#FFB984",
    "#3B9886",
    "#7A77AB",
    "#73B5E3",
    "#F28360",
    "#F7B846",
    "#97CD78",
    "#DCBFFC"
  )


arcadia_poppies <- list(
  color_dict = list(
    "arcadia:concord" = "#341E60",
    "arcadia:aegean" = "#5088C5",
    "arcadia:vitalblue" = "#73B5E3",
    "arcadia:paper" = "#FCFCFC",
    "arcadia:dress" = "#F8C5C1",
    "arcadia:amber" = "#F28360",
    "arcadia:dragon" = "#C85152",
    "arcadia:redwood" = "#52180A"
  ),
  values = c(0, 0.235, 0.33, 0.5, 0.68, 0.785, 1.0)
)

# Now, define some functions to make plotting using these colors as effective as possible:
# Function to select N maximally distinct colors from a palette
arcadia_color_discrete <- function(palette, n) {
  # Convert the palette to RGB color space
  rgb_palette <- grDevices::col2rgb(palette)

  # Convert RGB to HSV
  hsv_palette <- t(apply(rgb_palette, 2, function(color) {
    grDevices::rgb2hsv(matrix(color, nrow = 3))
  }))

  # Hue categories and their corresponding hue values
  hue_sequence <- c(
    "blue",
    "red-violet",
    "orange",
    "yellow-green",
    "blue-violet",
    "red",
    "yellow-orange",
    "green",
    "violet",
    "red-orange",
    "yellow",
    "blue-green"
  )
  hue_values <- c(
    240 / 360,
    300 / 360,
    30 / 360,
    90 / 360,
    270 / 360,
    0,
    45 / 360,
    120 / 360,
    270 / 360,
    15 / 360,
    60 / 360,
    210 / 360
  )

  # Categorize colors by hue
  categorize_by_hue <- function(h) {
    differences <- abs(hue_values - h)
    category_index <- which.min(differences)
    hue_sequence[category_index]
  }

  hue_categories <- sapply(hsv_palette[, 1], categorize_by_hue)

  # Sort colors within each hue category by saturation
  sorted_colors <- c()
  for (hue_cat in hue_sequence) {
    indices_in_category <- which(hue_categories == hue_cat)
    sorted_by_saturation <-
      indices_in_category[order(-hsv_palette[indices_in_category, 2])]
    sorted_colors <- c(sorted_colors, palette[sorted_by_saturation])
  }
  # Create an updated palette that includes the number of unique colors we need
  sorted_colors <- colorRampPalette(sorted_colors)(n)
  sorted_colors <- sorted_colors[-which(sorted_colors %in% palette)]

  # Loop through, pulling out colors that span the full color space,
  # removing these for the sorted list of colors, and repeating until we have a
  # new, resampled color vector that has sufficient separation between all
  # adjacent colors

  # Populate this with the starting palette - we want to use these colors first
  final_colors <- palette

  while (length(final_colors) < n) {
    # Figure out how many more colors we need to sample
    remaining_indices <- n - length(final_colors)

    # Sample at most 6 new colors that span the color space
    if (remaining_indices >= 6) {
      sampled_cols <-
        seq.int(
          from = 1,
          to = length(sorted_colors),
          length.out = 6
        )
    } else {
      sampled_cols <-
        seq.int(
          from = 1,
          to = length(sorted_colors),
          length.out = remaining_indices
        )
    }
    # Add to the list of final colors
    final_colors <- c(final_colors, sorted_colors[sampled_cols])
    # And remove from the vector of sorted colors
    sorted_colors <- sorted_colors[-sampled_cols]
  }

  return(final_colors)
}

# a function for min-max normalization that can
# handle the presence of NAs
minmax_norm <- function(x) {
  # Save the original positions of NA values
  na_positions <- is.na(x)

  # Ignore NA values while calculating min and max
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)

  # Normalize x, but ignore NA values
  x <- (x - min_x) / (max_x - min_x)

  # Replace the original positions of NA values
  x[na_positions] <- NA

  return(x)
}

# This function is called within the "run_phylo_profiling" function, and simply
# converts the dbscan kNN object (a list containing a matrix of distances),
# and a matrix of corresponding ids) into a dgCMatrix (compressed, sparse,
# and column-oriented). This is for efficiency and so that it may be used
# as input downstream.
convert_to_dgCMatrix <- function(idx, dist) {
  n <- nrow(idx) # Number of vertices
  n_neighbors <- ncol(idx) # Number of neighbors for each vertex

  # Prepare row indices, column indices, & distance values for the sparse matrix
  row_indices <- as.vector(t(idx))
  col_indices <- rep(1:n, each = n_neighbors)
  distance_values <- as.vector(t(dist))

  # Create the sparse matrix
  sparse_matrix <-
    sparseMatrix(
      i = row_indices,
      j = col_indices,
      x = distance_values,
      dims = c(n, n)
    )
  # Make the matrix symmetric
  sparse_matrix <- (sparse_matrix + t(sparse_matrix)) / 2

  return(sparse_matrix)
}
