# Run from analysis: Rscript "scripts/06_protected_area_overlap.R"
protected_area_overlap <- function(analysis = getwd(), release = "2026-10") {
  .libPaths(c(file.path(analysis, "scripts/.r-library"), .libPaths()))
  maxent <- file.path(analysis, "outputs/research question 5")
  output <- file.path(analysis, "outputs/research question 6")
  figures <- file.path(output, "figures")
  tables <- file.path(output, "tables")
  source_dir <- file.path(analysis, "Datasets/mapping/protected_planet")
  boundary <- file.path(analysis, "Datasets/mapping/boundaries/ken_admbnda_adm0_iebc_20191031.shp")
  source_files <- file.path(source_dir, paste0("WDPCA_", release, "_KEN_", c("polygons", "points"), ".geojson"))

  species <- c("Coffea fadenii", "Afrocarpus usambarensis")
  scenarios <- c("narrow_baseline", "broad_baseline", "shortlist")
  scenario_labels <- c("Narrow baseline", "Broad baseline", "Across-scenario shortlist")
  results <- expand.grid(scenario = scenarios, species = species, stringsAsFactors = FALSE)
  results <- results[c("species", "scenario")]
  rasters <- file.path(maxent, "figures", paste0(gsub(" ", "_", results$species),
    "_", results$scenario, "_q10_suitable.tif"))
  missing <- c(source_files, boundary, rasters)
  missing <- missing[!file.exists(missing)]
  if (length(missing)) stop("Missing inputs: ", paste(missing, collapse = "\n"))

  # Restrict the October 2026 WDPCA extract to protected areas in force.
  read_records <- function(path) {
    data <- sf::st_read(path, quiet = TRUE)
    stopifnot(all(c("site_type", "status", "site_id", "desig_eng", "desig_type") %in% names(data)))
    data[which(data$site_type == "PA" & !tolower(data$status) %in% c("proposed", "not reported")), ]
  }
  polygons <- read_records(source_files[1])
  points <- read_records(source_files[2])
  stopifnot(nrow(polygons) > 0, !any(sf::st_is_empty(polygons)),
    all(sf::st_geometry_type(polygons) %in% c("POLYGON", "MULTIPOLYGON")))
  point_only <- points[!points$site_id %in% polygons$site_id, ]

  # Equal-area geometry gives the fraction of each original cell inside the
  # dissolved Kenyan boundaries. Weight that fraction by Maxent's cellSize area.
  equal_area <- 6933
  kenya <- sf::st_union(sf::st_make_valid(sf::st_transform(sf::st_read(boundary, quiet = TRUE), equal_area)))
  protected <- sf::st_union(sf::st_make_valid(sf::st_transform(polygons, equal_area)))
  protected <- sf::st_intersection(protected, kenya)

  template <- terra::rast(rasters[1])
  masks <- lapply(rasters, terra::rast)
  reference <- terra::rast(file.path(maxent, "figures/Coffea_fadenii_reference_suitability.tif"))
  stopifnot(terra::compareGeom(template, reference), all(vapply(masks, function(mask)
    terra::compareGeom(template, mask), logical(1))))
  prediction_grid <- !is.na(terra::values(reference, mat = FALSE))
  values <- vapply(masks, terra::values, numeric(terra::ncell(template)), mat = FALSE)
  stopifnot(all(is.na(values) == !prediction_grid), all(values[!is.na(values)] %in% 0:1))
  suitable <- values == 1 & !is.na(values)
  cell_area <- terra::values(terra::cellSize(template, unit = "km"), mat = FALSE)

  cell_grid <- template
  ids <- seq_len(terra::ncell(template))
  terra::values(cell_grid) <- ifelse(rowSums(suitable) > 0, ids, NA_integer_)
  names(cell_grid) <- "cell"
  cells <- sf::st_transform(sf::st_as_sf(terra::as.polygons(cell_grid,
    aggregate = FALSE, values = TRUE, na.rm = TRUE)), equal_area)
  inside <- suppressWarnings(sf::st_intersection(cells, protected))
  fractions <- numeric(length(ids))
  if (nrow(inside)) {
    area_by_cell <- tapply(as.numeric(sf::st_area(inside)), inside$cell, sum)
    matched <- match(as.integer(names(area_by_cell)), cells$cell)
    fractions[cells$cell[matched]] <- area_by_cell / as.numeric(sf::st_area(cells[matched, ]))
  }
  stopifnot(all(is.finite(fractions)), all(fractions >= 0 & fractions <= 1 + 1e-7))
  results$suitable_km2 <- colSums(suitable * cell_area)
  results$suitable_in_protected_areas_km2 <- colSums(suitable * cell_area * pmin(fractions, 1))
  results$percent_protected <- 100 * results$suitable_in_protected_areas_km2 / results$suitable_km2

  baseline <- read.csv(file.path(maxent, "tables/scenario_suitable_area_km2.csv"))
  shortlist <- read.csv(file.path(maxent, "tables/prediction_coverage_and_agreement.csv"))
  expected <- vapply(seq_len(nrow(results)), function(i) {
    sp <- results$species[i]
    scenario <- results$scenario[i]
    if (scenario == "shortlist") return(shortlist$candidate_km2[shortlist$species == sp])
    baseline$q10_scenario_km2[baseline$species == sp & baseline$scenario == scenario]
  }, numeric(1))
  stopifnot(all(abs(results$suitable_km2 - expected) < .01))

  for (folder in c(figures, tables)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  # Format the manuscript table after calculating areas at full precision.
  table <- data.frame(
    Species = results$species,
    "Suitability layer" = scenario_labels[match(results$scenario, scenarios)],
    "Total area (km2)" = formatC(results$suitable_km2, format = "f", digits = 0, big.mark = ","),
    "Area inside protected areas (km2)" = formatC(results$suitable_in_protected_areas_km2,
      format = "f", digits = 0, big.mark = ","),
    "Inside (%)" = sprintf("%.1f", results$percent_protected),
    check.names = FALSE
  )
  csv <- capture.output(write.csv(table, row.names = FALSE))
  csv[1] <- gsub("km2", "km\u00b2", csv[1], fixed = TRUE, useBytes = TRUE)
  writeLines(csv, file.path(tables, "Table_S8_Protected_area_overlap_with_predicted_climatic_suitability.csv"), useBytes = TRUE)
  names(table) <- gsub("km2", "km\u00b2", names(table), fixed = TRUE, useBytes = TRUE)
  table$Species <- paste0("*", table$Species, "*")
  writeLines(c(paste0("| ", paste(names(table), collapse = " | "), " |"),
    "| --- | --- | ---: | ---: | ---: |",
    apply(table, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))),
    file.path(tables, "Table_S8_Protected_area_overlap_with_predicted_climatic_suitability.md"), useBytes = TRUE)
  write.csv(sf::st_drop_geometry(point_only), file.path(tables, "excluded_point_only_sites.csv"), row.names = FALSE)


  palette <- c(outside = "#D4B58A", inside = "#297A3F", boundary = "#332288")
  bounds <- sf::st_bbox(kenya)
  survey_sites <- read.csv(file.path(analysis,
    "Datasets/habitat_suitability/exploratory/forest_reference_points.csv"))
  taita <- sf::st_sfc(sf::st_point(c(mean(survey_sites$lon), mean(survey_sites$lat))), crs = 4326)
  taita_xy <- sf::st_coordinates(sf::st_transform(taita, equal_area))
  longitude <- seq(34, 42, by = 2)
  latitude <- seq(-4, 4, by = 2)
  longitude_x <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(data.frame(lon = longitude, lat = 0), coords = c("lon", "lat"), crs = 4326),
    equal_area))[, 1]
  latitude_y <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(data.frame(lon = 38, lat = latitude), coords = c("lon", "lat"), crs = 4326),
    equal_area))[, 2]
  draw_maps <- function() {
    layout(rbind(c(1, 0, 2), c(3, 3, 3)),
      widths = c(1, lcm(.2), 1), heights = c(1, lcm(1.28)))
    par(bg = "white", col = "black", col.axis = "black")
    for (i in which(results$scenario == "shortlist")) {
      panel <- match(results$species[i], species)
      par(mai = c(.23, .32, .05, .06))
      map_height <- bounds[4] - bounds[2]
      plot(kenya, col = "#FCFCFC", border = NA, xlim = bounds[c(1, 3)],
        ylim = bounds[c(2, 4)] + c(-.025, .12) * map_height, yaxs = "i")
      plot(sf::st_geometry(cells[suitable[cells$cell, i], ]), col = palette["outside"], border = NA, add = TRUE)
      # These are the same intersections used to calculate fractions above.
      if (nrow(inside)) plot(sf::st_geometry(inside[suitable[inside$cell, i], ]),
        col = palette["inside"], border = NA, add = TRUE)
      plot(protected, border = palette["boundary"], lwd = .4, add = TRUE)
      plot(kenya, border = "#606060", lwd = .85, add = TRUE)
      points(taita_xy, pch = 21, bg = "white", col = "black", cex = .95)
      text(taita_xy, labels = "Taita Hills", pos = 4, offset = .40, cex = .75 + 1/12)
      limits <- par("usr")
      heading_x <- mean(limits[1:2])
      text(heading_x, bounds[4] + .096 * map_height,
        bquote(bold(.(paste0("(", letters[panel], ")"))) ~ italic(.(results$species[i]))), cex = 1 + 2/12)
      text(heading_x, bounds[4] + .05 * map_height,
        sprintf("%s km\u00b2 inside protected areas (%.1f%%)",
          format(round(results$suitable_in_protected_areas_km2[i]), big.mark = ",", trim = TRUE),
          results$percent_protected[i]), cex = .77 + 4/12)
      compass_x <- limits[1] + .93 * diff(limits[1:2])
      compass_y <- limits[3] + .5 * diff(limits[3:4])
      compass_radius <- .047 * diff(limits[1:2]) + diff(limits[1:2]) / (144 * par("pin")[1])
      arc <- seq(.62 * pi, 2.38 * pi, length.out = 120)
      lines(compass_x + compass_radius * cos(arc),
        compass_y + compass_radius * sin(arc), lwd = 1.25)
      needle_y <- compass_y + compass_radius * c(.67, -.45, -.22)
      polygon(compass_x + compass_radius * c(0, -.43, 0), needle_y,
        col = "white", border = "black", lwd = 1.1)
      polygon(compass_x + compass_radius * c(0, .43, 0), needle_y,
        col = "black", border = "black", lwd = 1.1)
      text(compass_x, compass_y + compass_radius, "N", font = 2, cex = .78 + 1/12)
      # EPSG:6933 horizontal scale at the equator.
      scale_right <- limits[1] + .94 * diff(limits[1:2])
      scale_left <- scale_right - 100000 * cos(pi / 6)
      scale_y <- limits[3] + .06 * diff(limits[3:4])
      scale_mid <- mean(c(scale_left, scale_right))
      bar_height <- (.008 + 1 / (72 * par("pin")[2])) * diff(limits[3:4])
      rect(c(scale_left, scale_mid), scale_y, c(scale_mid, scale_right),
        scale_y + bar_height, col = c("black", "white"), border = "black")
      text(scale_mid, scale_y + 1.5 * bar_height, "100 km",
        adj = c(.5, 0), cex = .78, font = 2)
      box(col = "grey50")
      axis(1, at = longitude_x, labels = paste0(longitude, "\u00b0 E"),
        cex.axis = .72 + 1/12, mgp = c(1, .35, 0), tcl = -.22, col = "grey50", col.ticks = "grey50")
      axis(2, at = latitude_y, labels = paste0(abs(latitude), "\u00b0",
        ifelse(latitude < 0, " S", ifelse(latitude > 0, " N", ""))),
        las = 1, cex.axis = .72 + 1/12, mgp = c(1, .35, 0), tcl = -.22, col = "grey50", col.ticks = "grey50")
    }
    par(mar = rep(0, 4))
    plot.new()
    plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    shortlist_key <- legend(.5, .99, c("Shortlisted area outside protected areas",
      "Shortlisted area within protected areas"), fill = palette[1:2],
      border = NA, horiz = TRUE, xjust = .5, bty = "n", text.col = "black", cex = .95 + 2/12)
    legend(shortlist_key$rect$left, .54, "Mapped protected-area boundary",
      lty = 1, lwd = .8, col = palette["boundary"], seg.len = 2.2,
      xjust = 0, bty = "n", text.col = "black", cex = .95 + 2/12)
    second_key_left <- shortlist_key$rect$left + diff(shortlist_key$text$x)
    legend(second_key_left, .54, "Taita Hills", pch = 21, pt.bg = "white",
      pt.cex = 1 + 2/12, xjust = 0, bty = "n", text.col = "black", cex = .95 + 2/12)
  }

  width <- 6.25 * 1.15
  height <- 5.18
  path <- file.path(figures, "Figure_8_Protected_area_overlap_with_the_climatic_suitability_shortlist")
  grDevices::cairo_pdf(paste0(path, ".pdf"), width = width, height = height, family = "Arial")
  tryCatch(draw_maps(), finally = grDevices::dev.off())
  grDevices::png(paste0(path, ".png"), width = round(width * 300), height = height * 300,
    res = 300, type = "cairo", family = "Arial")
  tryCatch(draw_maps(), finally = grDevices::dev.off())
  invisible(results)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  analysis <- if (length(script) && file.exists(script[1])) {
    dirname(dirname(normalizePath(script[1], mustWork = TRUE)))
  } else {
    working_dir <- getwd()
    if (file.exists(file.path(working_dir, "scripts", "06_protected_area_overlap.R"))) {
      working_dir
    } else if (file.exists(file.path(working_dir, "analysis", "scripts", "06_protected_area_overlap.R"))) {
      file.path(working_dir, "analysis")
    } else {
      stop("Set the working directory to the project root or its analysis folder, or call protected_area_overlap() with the analysis folder path.")
    }
  }
  print(protected_area_overlap(analysis, release = if (length(args)) args[1] else "2026-10"))
  source(file.path(analysis, 'scripts', 'combine_output_tables.R'), local = TRUE)
  combine_output_tables(file.path(analysis, 'outputs', 'research question 6'), '06_protected_area_overlap.R')
}
