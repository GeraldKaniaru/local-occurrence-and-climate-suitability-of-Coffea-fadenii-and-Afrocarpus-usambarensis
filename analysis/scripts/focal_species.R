# ============================================================
# 0. Focal species
# ============================================================

# Run with Rscript analysis/scripts/focal_species.R from the project.
focal_species <- function(root) {
  if (!requireNamespace("jpeg", quietly = TRUE)) stop("Install the jpeg package to read the focal-species photos.")
  # Coffea: foliage, fruit, plant. Afrocarpus: foliage, bark, crown.
  sources <- file.path(root, "analysis", "focal_species_photos", c(
    "20210623_121107.jpg", "20210623_121043.jpg", "20210622_103625.jpg",
    "20210623_104226.jpg", "20210623_104500.jpg", "20210624_105033.jpg"
  ))
  missing <- sources[!file.exists(sources)]
  if (length(missing)) stop("Missing focal-species photos:\n", paste(missing, collapse = "\n"))
  photos <- lapply(sources, jpeg::readJPEG)
  tile_widths <- 5 * c(205, 252, 167, 215, 245, 164)
  row_heights <- 5 * c(217, 227)
  divider <- 10
  width <- sum(tile_widths[1:3])
  height <- sum(row_heights) + divider
  out <- file.path(root, "analysis", "outputs", "focal_species", "figures")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(out, "Figure_S1_Focal_species_in_the_Taita_Hills.png")
  grDevices::png(target, width = width, height = height, res = 300)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "#BDBDBD", col = NA))
  for (i in seq_along(photos)) {
      row <- if (i <= 3L) 1L else 2L
      first <- if (row == 1L) 1L else 4L
      left <- if (i == first) 0 else sum(tile_widths[first:(i - 1L)])
      bottom <- if (row == 1L) row_heights[2] + divider else 0
      tile_width <- tile_widths[i]
      tile_height <- row_heights[row]
      photo <- photos[[i]]
      # Centre-crop each photo to fill its panel while preserving aspect ratio.
      scale <- max(tile_width / ncol(photo), tile_height / nrow(photo))
      crop_width <- min(ncol(photo), round(tile_width / scale))
      crop_height <- min(nrow(photo), round(tile_height / scale))
      x_start <- (ncol(photo) - crop_width) %/% 2L
      y_start <- (nrow(photo) - crop_height) %/% 2L
      photo <- photo[y_start + seq_len(crop_height), x_start + seq_len(crop_width), , drop = FALSE]
      grid::grid.raster(photo, x = (left + tile_width / 2) / width,
        y = (bottom + tile_height / 2) / height,
        width = tile_width / width, height = tile_height / height)
      x <- (left + 30) / width
      y <- (bottom + tile_height - 28) / height
      grid::grid.rect(x, y, width = grid::unit(9, "mm"), height = grid::unit(7, "mm"),
        just = c("left", "top"), gp = grid::gpar(fill = "white", col = NA))
      grid::grid.text(paste0("(", letters[i], ")"),
        x = grid::unit(x, "npc") + grid::unit(4.5, "mm"),
        y = grid::unit(y, "npc") - grid::unit(3.5, "mm"),
        gp = grid::gpar(fontface = "bold", fontsize = 13))
  }
  message("Photo plate saved to ", target)
  invisible(target)
}

# ============================================================
# 1. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  root <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
  while (!file.exists(file.path(root, "analysis/scripts/focal_species.R"))) {
    parent <- dirname(root)
    if (identical(root, parent)) stop("Cannot locate the Taita project.")
    root <- parent
  }
  focal_species(root)
}
