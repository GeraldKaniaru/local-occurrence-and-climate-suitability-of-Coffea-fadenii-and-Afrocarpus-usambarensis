# ============================================================
# 0. Setup
# ============================================================

# Reproducible study-area plate; local inputs only, source coordinates unchanged.
# Rscript -- analysis/scripts/study_area.R

# ============================================================
# 1. Run study area
# ============================================================

run_study_area <- function(root) {
  .libPaths(c(file.path(root, 'analysis/scripts/.r-library'), .libPaths()))
  stopifnot(requireNamespace('terra', quietly = TRUE), requireNamespace('readxl', quietly = TRUE))
  out <- file.path(root, 'analysis/outputs/study_area')
  figures <- file.path(out, 'figures')
  tables <- file.path(out, 'tables')
  for (path in c(figures, tables))dir.create(path, recursive = TRUE, showWarnings = FALSE)
  read <- function(p) read.csv(file.path(root, p), stringsAsFactors = FALSE)
  derived <- file.path(root, 'analysis/Datasets/derived')
  dir.create(derived, recursive = TRUE, showWarnings = FALSE)
  location_inputs <- 'analysis/Datasets/derived/worksheet_locations.csv'
  # Audit original workbook headers; agreement confirms transcription, not GPS accuracy.
  raw <- list();
  k <- 0
  for (f in c('Ngangao', 'Mbololo')) {
    workbook <- file.path(root, paste0('analysis/Datasets/field/', f, '.Mapping.Data.Franklinia.xlsx'))
    for (sh in readxl::excel_sheets(workbook)) {
      cells <- as.character(unlist(readxl::read_excel(workbook, sheet = sh, range = 'A1:H6', col_names = FALSE, .name_repair = 'minimal')))
      grab <- function(word) {v <- cells[!is.na(cells) & grepl(word, cells, ignore.case = TRUE)];
        if (length(v)) v[1] else NA_character_}
      parse_dm <- function(x, sign) {if (is.na(x)) return(NA_real_);
        z <- regmatches(x, gregexpr('[0-9]+(?:\\.[0-9]+)?', x, perl = TRUE))[[1]];
      if (length(z)<2) return(NA_real_);
      sign * (as.numeric(z[1]) + as.numeric(z[2]) / 60)}
      lat <- grab('Latitude');
      lon <- grab('Longitude');
      k <- k + 1
      raw[[k]] <- data.frame(plot_id = sprintf('%s_%02d', f, as.integer(sub('Plot.', '', sh, fixed = TRUE))), source_workbook = substring(workbook, nchar(root) + 2), sheet = sh, raw_latitude = lat, raw_longitude = lon, source_lat = parse_dm(lat, - 1), source_lon = parse_dm(lon, 1))
    }
  }
  worksheet_coordinates <- do.call(rbind, raw)
  write.csv(worksheet_coordinates, file.path(tables, 'worksheet_coordinate_transcription.csv'), row.names = FALSE)
  write.csv(worksheet_coordinates, file.path(derived, 'worksheet_coordinate_transcription.csv'), row.names = FALSE)
  available <- file.exists(file.path(root, location_inputs))
  env <- file.path(root, 'analysis/Datasets/mapping/environment')
  raster_inputs <- file.path(env, c('elevation.tif', 'tree_cover.tif'))
  missing <- raster_inputs[!file.exists(raster_inputs)]
  if (!any(available))missing <- c(paste('Survey ledger:', file.path(root, location_inputs)), missing)
  if (length(missing)) {
    png_name <- 'Figure_1_Study_areas_in_the_Taita_Hills.png'
    png_path <- file.path(figures, png_name)
    saved_png <- file.path(root, 'submission/Figures', png_name)
    restored <- !isTRUE(file.size(png_path)>0) && isTRUE(file.size(saved_png)>0) &&
    file.copy(saved_png, png_path, overwrite = TRUE)
    saved_figures <- list.files(figures, pattern = '^Figure_1_Study_areas_in_the_Taita_Hills(_render)?\\.(png|pdf)$', full.names = TRUE)
    figure_status <- if (restored) {
      'The saved manuscript figure was copied to the output folder.'
    } else if (any(file.size(saved_figures)>0, na.rm = TRUE)) {
      'Existing figure files were retained.'
    } else {
      'No existing figure files were found.'
    }
    status <- c('Study-area map not regenerated.', figure_status,
      sprintf('Transcribed coordinates from %d worksheets (%d with both coordinates).',
        nrow(worksheet_coordinates), sum(is.finite(worksheet_coordinates$source_lat) & is.finite(worksheet_coordinates$source_lon))),
      'Worksheet transcriptions do not establish eligibility or cultivated status.',
      'Missing inputs:', missing)
    message(paste(status, collapse = '\n'))
    return(invisible(FALSE))
  }
  q <- read(location_inputs[which(available)[1]])
  r <- q[q$primary_eligible & !q$cultivated & is.finite(q$lon) & is.finite(q$lat),]
  stopifnot(nrow(r) == 23, all(abs(r$lat)<90), all(abs(r$lon)<180))
  audit <- merge(q, worksheet_coordinates, by = 'plot_id', all.x = TRUE, sort = FALSE)
  audit$transcription_matches <- abs(audit$lat - audit$source_lat)<1e-8 & abs(audit$lon - audit$source_lon)<1e-8
  stopifnot(all(audit$transcription_matches))
  audit$verification_priority <- ifelse(audit$plot_id == 'Ngangao_06', 'Ngangao assignment confirmed by user; coordinate unchanged; southern position versus northern-side source description remains a coordinate audit issue', 'GPS accuracy, datum and forest assignment require field confirmation')
  write.csv(audit, file.path(tables, 'coordinate_source_audit.csv'), row.names = FALSE)
  # Retain site-reference provenance but do not use those independent points on the map.
  sites <- read('analysis/Datasets/mapping/Study sites.csv');
  sites <- sites[match(c('Ngangao forest', 'Mbololo forest'), sites$study_site),];
  stopifnot(!anyNA(sites$study_site));
  sites$map_use <- 'Panel c only: approximate forest reference points, not worksheet records'
  write.csv(sites, file.path(tables, 'study_site_coordinate_sources.csv'), row.names = FALSE)
  kenya <- terra::vect(file.path(root, 'analysis/Datasets/mapping/boundaries/ken_admbnda_adm0_iebc_20191031.shp'))
  counties <- terra::vect(file.path(root, 'analysis/Datasets/mapping/boundaries/ken_admbnda_adm1_iebc_20191031.shp'))
  taita <- counties[counties$ADM1_EN == 'Taita Taveta',]
  africa <- terra::project(terra::vect(file.path(root, 'analysis/Datasets/mapping/boundaries/Africa_Boundaries.shp')), 'EPSG:4326')
  # Fixed, shared display extent; all records included with margins. No smoothing.
  detail <- c(38.30, 38.50, - 3.42, - 3.29)
  terrain_extent <- c(38.23, 38.54, - 3.46, - 3.25)
  stopifnot(all(r$lon>detail[1] & r$lon<detail[2] & r$lat>detail[3] & r$lat<detail[4]))
  elev <- terra::crop(terra::rast(file.path(env, 'elevation.tif')), terra::ext(terrain_extent), snap = 'out')
  trees <- 100 * terra::crop(terra::rast(file.path(env, 'tree_cover.tif')), terra::ext(detail), snap = 'out')
  stopifnot(max(abs(terra::res(elev) - 1 / 120))<1e-8, max(abs(terra::res(trees) - 1 / 120))<1e-8)
  r$environmental_cell <- terra::cellFromXY(trees, as.matrix(r[, c('lon', 'lat')]))
  r$records_in_same_cell <- as.integer(table(r$environmental_cell)[as.character(r$environmental_cell)])
  write.csv(r, file.path(tables, 'mapped_primary_survey_records.csv'), row.names = FALSE)
  # Geodesic diagnostics are checks, not corrections or boundary verification.
  distances <- as.matrix(terra::distance(as.matrix(r[, c('lon', 'lat')]), lonlat = TRUE))
  diag(distances) <- Inf
  r$nearest_other_record_m <- apply(distances, 1, min)
  write.csv(r, file.path(tables, 'mapped_primary_survey_records.csv'), row.names = FALSE)
  refs <- sites
  refs$nearest_same_forest_record_m <- vapply(seq_len(nrow(refs)), function(i) {
      f <- sub(' forest', '', refs$study_site[i]);
      z <- r[r$forest == f,]
      min(terra::distance(matrix(c(refs$long[i], refs$lat[i]), ncol = 2), as.matrix(z[, c('lon', 'lat')]), lonlat = TRUE))
    }, numeric(1))
  write.csv(refs, file.path(tables, 'study_site_coordinate_sources.csv'), row.names = FALSE)
  # Precompute polygon rings once; inset simplification affects display only.
  rings <- function(v, tolerance = 0) {if (tolerance>0)v <- terra::simplifyGeom(v, tolerance = tolerance);
    m <- terra::geom(v);
    split(data.frame(x = m[, 'x'], y = m[, 'y']), interaction(m[, 'geom'], m[, 'part'], drop = TRUE))}
  kg <- rings(kenya);
  cg <- rings(counties);
  tg <- rings(taita);
  ag <- rings(africa, .03)
  polys <- function(g, fill, border = 'grey40', lwd = .6)for (x in g)polygon(x$x, x$y, col = fill, border = border, lwd = lwd)
  cultivated <- q[q$cultivated & is.finite(q$lon) & is.finite(q$lat),]
  stopifnot(nrow(cultivated) == 2, all(cultivated$lon>detail[1] & cultivated$lon<detail[2] & cultivated$lat>detail[3] & cultivated$lat<detail[4]))
  display_records <- rbind(q[match(r$plot_id, q$plot_id),], cultivated)
  display_records$display_context <- ifelse(display_records$cultivated, 'Cultivated', 'Primary forest-associated')
  write.csv(display_records, file.path(tables, 'mapped_display_records.csv'), row.names = FALSE)
  county_extent <- as.vector(terra::ext(taita)) + c( - .035, .035, - .035, .035)
  colours <- c(Ngangao = '#D4B58A', Mbololo = '#508A92')
  # Give each map an exact extent and a plot rectangle with geographic aspect.
  frame <- function(fig, extent, new = TRUE, margins = c(.55, .65, .35, .8)) {
    par(fig = fig, new = new, mar = rep(0, 4), family = 'Arial', las = 1)
    fin <- par('fin');
    avail <- fin - c(sum(margins[c(2, 4)]), sum(margins[c(1, 3)]))
    ratio <- diff(extent[1:2]) * cos(mean(extent[3:4]) * pi / 180) / diff(extent[3:4])
    w <- min(avail[1], avail[2] * ratio);
    h <- w / ratio
    left <- (margins[2] + (avail[1] - w) / 2) / fin[1];
    bottom <- (margins[1] + (avail[2] - h) / 2) / fin[2]
    par(plt = c(left, left + w / fin[1], bottom, bottom + h / fin[2]))
    plot.new();
    plot.window(xlim = extent[1:2], ylim = extent[3:4], xaxs = 'i', yaxs = 'i')
    rect(extent[1], extent[3], extent[2], extent[4], col = 'white', border = NA)
  }
  # Geodesic scale at its drawn latitude; verify endpoints with terra's geodesic distance.
  scale_checks <- list()
  scale_bar <- function(km, x, y) {
    d <- terra::distance(matrix(c(x, y), ncol = 2), matrix(c(x + 1, y), ncol = 2), lonlat = TRUE)[1]
    len <- km * 1000 / d;
    h <- diff(par('usr')[3:4]) * .014
    rect(x, y, x + len / 2, y + h, col = 'black');
    rect(x + len / 2, y, x + len, y + h, col = 'white')
    measured <- terra::distance(matrix(c(x, y), ncol = 2), matrix(c(x + len, y), ncol = 2), lonlat = TRUE)[1]
    stopifnot(abs(measured / (km * 1000) - 1)<.001)
    scale_checks[[length(scale_checks) + 1]] <<- data.frame(requested_m = km * 1000, geodesic_m = measured, latitude = y)
    c(x + len / 2, y + 3 * h)
  }
  raster_map <- function(z, palette, limits) {
    m <- as.matrix(z, wide = TRUE);
    u <- par('usr');
    rect(u[1], u[3], u[2], u[4], col = '#dadada', border = NA)
    image(terra::xFromCol(z, 1:ncol(z)), rev(terra::yFromRow(z, 1:nrow(z))), t(m[nrow(m):1,, drop = FALSE]), add = TRUE, col = palette, zlim = limits)
    # Compact colour bar; endpoints explicitly labelled.
    dx <- diff(u[1:2]);
    dy <- diff(u[3:4]);
    x <- u[2] + .035 * dx;
    w <- .022 * dx;
    yy <- seq(u[3], u[4], length.out = length(palette) + 1)
    rect(x, head(yy, - 1), x + w, tail(yy, - 1), col = palette, border = NA, xpd = NA)
    rect(x, u[3], x + w, u[4], xpd = NA)
    invisible(list(right = x + w, dx = dx, bottom = u[3], top = u[4], dy = dy, center = mean(u[1:2])))
  }
  records_layer <- function() {
    for (f in names(colours)) {z <- r[r$forest == f,];
      points(z$lon, z$lat, pch = 21, bg = colours[f], col = 'black', lwd = .8)}
    points(cultivated$lon, cultivated$lat, pch = 24, bg = '#CC79A7', col = 'black', lwd = .8)
  }
  plate <- function() {
    par(family = 'Arial', bg = 'white', col = 'black')
    # Panel (a) Kenya and Africa inset
    e <- c(33.5, 42.4, - 5.1, 5.9)
    frame(c(0, .39, .46, 1), e, FALSE, margins = c(.5, .65, .15, .12))
    polys(kg, '#f6f6f3', 'black', .8);
    polys(cg, NA, '#cccccc', .4);
    polys(tg, '#bdd5b1', '#416b47', 1.1)
    be <- county_extent;
    rect(be[1], be[3], be[2], be[4], border = '#6b4c9a', lwd = 1.4, lty = 2)
    text(38.1, - 2.3, 'Taita-Taveta', cex = 1.1, font = 2, col = 'black')
    text(39.5, - 4.2, 'b', cex = 1.1, font = 2, col = 'black')
    scale_label <- scale_bar(200, 34.2, - 4.65)
    text(scale_label[1], scale_label[2], '200 km', cex = 1.1, font = 2, col = 'black')
    axis(1, at = c(34, 38, 42), labels = c('34\u00b0 E', '38\u00b0 E', '42\u00b0 E'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    axis(2, at = c( - 4, 0, 4), labels = c('4\u00b0 S', '0\u00b0 N', '4\u00b0 N'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    box(col = '#999999')
    panel_limits <- par('usr')
    text(panel_limits[2] - .025 * diff(panel_limits[1:2]),
      panel_limits[4] - .035 * diff(panel_limits[3:4]),
      '(a)', adj = c(1, 1), cex = 1.9, font = 2, col = 'black')
    # North arrow at the right-center of the Kenya panel.
    arrows(41.8, - .35, 41.8, .65, length = .1, lwd = 1.4, col = 'black')
    text(41.8, .92, 'N', cex = 1.1, font = 2, col = 'black')
    # Locate inset using Kenya plotting-region device coordinates, keeping it inside.
    fx <- grconvertX(c(33.7, 36.6), from = 'user', to = 'ndc');
    fy <- grconvertY(c(2.65, 5.85), from = 'user', to = 'ndc')
    frame(c(fx, fy), c( - 20, 55, - 36, 38), TRUE, margins = rep(.015, 4))
    polys(ag, 'white', '#555555', .4);
    polys(kg, '#29854a', 'black', .6);
    box(col = '#aaaaaa')
    # Panel (b) Taita-Taveta County
    frame(c(0, .39, .07, .46), be, margins = c(.5, .65, .12, .12))
    polys(tg, '#d5e4ca', '#555555', .8)
    rect(terrain_extent[1], terrain_extent[3], terrain_extent[2], terrain_extent[4], border = '#6b4c9a', lwd = 1.5, lty = 2)
    text(38.5, - 3.15, 'c', col = 'black', cex = 1.1, font = 2)
    text(38.05, - 2.85, 'Taita-Taveta', cex = 1.1, font = 2, col = 'black')
    text(38.25, - 3.65, 'Taita Hills', cex = 1.1, font = 2, col = 'black')
    scale_label <- scale_bar(30, 37.70, - 4.10)
    text(scale_label[1], scale_label[2], '30 km', cex = 1.1, font = 2, col = 'black')
    axis(1, at = c(37.7, 38.1, 38.5, 39.0), labels = c('37.7\u00b0 E', '38.1\u00b0 E', '38.5\u00b0 E', '39.0\u00b0 E'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    axis(2, at = c( - 4, - 3.5, - 3), labels = c('4.0\u00b0 S', '3.5\u00b0 S', '3.0\u00b0 S'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    box(col = '#999999')
    panel_limits <- par('usr')
    text(panel_limits[2] - .025 * diff(panel_limits[1:2]),
      panel_limits[4] - .035 * diff(panel_limits[3:4]),
      '(b)', adj = c(1, 1), cex = 1.9, font = 2, col = 'black')
    # Panel (c) Elevation and forest reference locations
    frame(c(.395, 1, .525, 1), terrain_extent, margins = c(.5, .85, .6, 1.0))
    bar <- raster_map(elev, hcl.colors(55, 'Terrain'), c(400, 2200))
    ticks <- sort(unique(c(400, 2200, pretty(c(400, 2200), 4))))
    ticks <- ticks[ticks >= 400 & ticks <= 2200]
    text(bar$right + .015 * bar$dx, bar$bottom + (ticks - 400) / (2200 - 400) * bar$dy,
      ticks, adj = 0, cex = 1.2, font = 1, col = 'black', xpd = NA)
    text(bar$center, bar$top + .065 * bar$dy, 'Elevation (m)',
      cex = 1.25, font = 2, col = 'black', xpd = NA)
    terra::contour(elev, add = TRUE, levels = seq(400, 2200, 400), col = 'black', labcex = 1.15, lwd = .65)
    points(sites$long, sites$lat, pch = 21, bg = '#f7d65e', col = 'black')
    text(sites$long[1] - .014, sites$lat[1] - .016, 'Ngangao', col = 'black', cex = 1.1, font = 2)
    text(sites$long[2] + .006, sites$lat[2] + .018, 'Mbololo', col = 'black', cex = 1.1, font = 2)
    rect(detail[1], detail[3], detail[2], detail[4], border = '#6b4c9a', lty = 2, lwd = 1.3)
    text(detail[1] + .005, detail[4] - .005, 'd', adj = c(0, 1), cex = 1.1, font = 2, col = 'black')
    arrows(38.52, - 3.315, 38.52, - 3.295, length = .1, lwd = 1)
    text(38.52, - 3.288, 'N', cex = 1.1, font = 2, col = 'black')
    legend('bottomright', 'Approximate forest reference', pch = 21, pt.bg = '#f7d65e', pt.lwd = 1.5, bty = 'n', text.col = 'black', cex = 1.1, text.font = 2, inset = c(.035, .025), xpd = FALSE)
    scale_label <- scale_bar(5, 38.245, - 3.45)
    text(scale_label[1], scale_label[2], '5 km', cex = 1.1, font = 2, col = 'black')
    axis(1, at = c(38.25, 38.35, 38.45, 38.50), labels = c('38.25\u00b0 E', '38.35\u00b0 E', '38.45\u00b0 E', '38.50\u00b0 E'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    axis(2, at = c( - 3.45, - 3.40, - 3.35, - 3.30), labels = c('3.45\u00b0 S', '3.40\u00b0 S', '3.35\u00b0 S', '3.30\u00b0 S'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    box(col = '#999999')
    panel_limits <- par('usr')
    text(panel_limits[2] - .025 * diff(panel_limits[1:2]),
      panel_limits[4] - .035 * diff(panel_limits[3:4]),
      '(c)', adj = c(1, 1), cex = 1.9, font = 2, col = 'black')
    # Panel (d) Tree class fraction and survey locations
    frame(c(.395, 1, .085, .525), detail, margins = c(.5, .85, .35, 1.0))
    bar <- raster_map(trees, colorRampPalette(c('#f6f2df', '#b9d5a3', '#23643c'))(60), c(0, 100))
    ticks <- sort(unique(c(0, 100, pretty(c(0, 100), 4))))
    ticks <- ticks[ticks >= 0 & ticks <= 100]
    text(bar$right + .015 * bar$dx, bar$bottom + (ticks - 0) / (100 - 0) * bar$dy,
      ticks, adj = 0, cex = 1.2, font = 1, col = 'black', xpd = NA)
    text(bar$center, bar$top + .065 * bar$dy, 'Tree class fraction (%)',
      cex = 1.25, font = 2, col = 'black', xpd = NA)
    records_layer()
    scale_label <- scale_bar(3, 38.308, - 3.414)
    text(scale_label[1], scale_label[2], '3 km', cex = 1.1, font = 2, col = 'black')
    axis(1, at = c(38.30, 38.35, 38.40, 38.45, 38.50), labels = c('38.30\u00b0 E', '38.35\u00b0 E', '38.40\u00b0 E', '38.45\u00b0 E', '38.50\u00b0 E'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    axis(2, at = c( - 3.42, - 3.38, - 3.34, - 3.30), labels = c('3.42\u00b0 S', '3.38\u00b0 S', '3.34\u00b0 S', '3.30\u00b0 S'), cex.axis = 1.3, font.axis = 1, col.axis = 'black', tcl = - .25)
    box(col = '#999999')
    panel_limits <- par('usr')
    text(panel_limits[2] - .025 * diff(panel_limits[1:2]),
      panel_limits[4] - .035 * diff(panel_limits[3:4]),
      '(d)', adj = c(1, 1), cex = 1.9, font = 2, col = 'black')
    # Shared survey legend
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4));
    plot.new()
    legend(.70, .050, legend = c('Ngangao (n = 9)', 'Mbololo (n = 14)', 'Cultivated (n = 2)'),
      pch = c(21, 21, 24), pt.bg = c(colours, '#CC79A7'), col = c('white', 'white', 'black'),
      bty = 'n', cex = 1.5, text.col = 'black', horiz = TRUE, xjust = .5, yjust = .5, x.intersp = .6)
  }
  stem <- file.path(figures, 'Figure_1_Study_areas_in_the_Taita_Hills')
  png(paste0(stem, '_render.png'), width = 4200, height = 3600, res = 300, pointsize = 12, type = 'cairo', family = 'Arial');
  plate();
  dev.off()
  cairo_pdf(paste0(stem, '_render.pdf'), width = 14, height = 12, pointsize = 12, family = 'Arial');
  plate();
  dev.off()
  for (ext in c('.png', '.pdf')) {
    target <- paste0(stem, ext);
    rendered <- paste0(stem, '_render', ext)
    if (file.exists(target))stopifnot(unlink(target) == 0L)
    stopifnot(file.rename(rendered, target))
  }
  write.csv(do.call(rbind, scale_checks), file.path(tables, 'scale_bar_verification.csv'), row.names = FALSE)
  write.csv(data.frame(panel = c('a', 'b', 'c', 'd'), xmin = c(33.5, county_extent[1], terrain_extent[1], detail[1]), xmax = c(42.4, county_extent[2], terrain_extent[2], detail[2]), ymin = c( - 5.1, county_extent[3], terrain_extent[3], detail[3]), ymax = c(5.9, county_extent[4], terrain_extent[4], detail[4])), file.path(tables, 'panel_extents.csv'), row.names = FALSE)
  metadata <- data.frame(layer = c('Elevation', 'Tree-class fraction'), source = c('WorldClim 2.1; SRTM-derived elevation', 'geodata WorldCover_trees_30s; documentation cites ESA WorldCover 2020 v100'), native_source_resolution = c('SRTM source; supplied WorldClim grid 30 arc-seconds', 'ESA WorldCover 10 m / 0.3 arc-seconds'), supplied_and_displayed_resolution = '30 arc-seconds (~925 m at these latitudes)', local_processing = c('Crop only; no smoothing or resampling', 'Crop and multiply provider fraction by 100; no smoothing, aggregation or resampling'), missing_cells_in_crop = c(sum(is.na(terra::values(elev))), sum(is.na(terra::values(trees)))), missing_display = 'Grey; no interpolation or imputation', method_limit = c('Coarse topography only', 'Provider fractions of cell classified as trees; fine-pixel weighting and missing-data denominator not specified in retrieved documentation'), documentation = c('https://worldclim.org/data/worldclim21.html', 'https://rspatial.github.io/geodata/reference/landcover.html'))
  write.csv(metadata, file.path(tables, 'environmental_layer_metadata.csv'), row.names = FALSE)
  caption <- paste('Study area and provisional survey-record locations. (a) Kenya, with Taita-Taveta County shaded and an Africa inset highlighting Kenya. The purple dashed rectangle marks the exact geographic extent of panel (b). (b) Taita-Taveta County and Taita Hills; the purple dashed rectangle marks the extent of panel (c). (c) Restored coloured WorldClim 2.1 SRTM-derived elevation overview, in metres, with coarse 400-m contours. Yellow circles are approximate forest reference locations from analysis/Datasets/mapping/Study sites.csv, not worksheet records or forest centroids. The purple dashed rectangle marks the exact extent of panel (d). (d) Percentage of each 30-arc-second cell classified as the tree-cover class in the geodata WorldCover-derived fraction raster; provider documentation cites ESA WorldCover 2020 v100, whose source classification has 10-m resolution. The provider aggregation represents class fraction, not canopy percentage; exact fine-pixel weighting and missing-pixel rules are not specified in the retrieved documentation. Both environmental layers are displayed on their supplied 30-arc-second grids (approximately 925 m), without local aggregation, resampling or smoothing. Missing cells would be grey; no values are imputed, and there are no missing cells in the displayed crops. The terrain colour scale progresses from green at lower elevations through yellow and peach to pale high-elevation cells; white or pale colours do not indicate missing data. Panel (d) plots the 23 primary forest-associated worksheet coordinates and two cultivated Mbololo records unchanged and without jitter. Tan circles mark Ngangao and teal circles mark Mbololo, pink triangles identify cultivated records. Ngangao_06 is shown with the same tan symbol as the other Ngangao records, consistent with its forest assignment and user confirmation. The cultivated records remain excluded from forest analyses. Symbols may overlap. Panel (c) reference points have a separate source and differ from the survey-record coordinates in panel (d); neither source has been independently field-verified. The axes use geographic WGS84 (EPSG:4326); this datum is confirmed for the GIS layers but only a computational assumption for the field records, whose original datum and accuracy remain unverified. Scale bars are checked against WGS84 geodesic distances at their drawn latitudes. Environmental cells provide landscape context and establish neither indigenous forest identity nor plot-level habitat conditions. Forest-fragment boundaries, plot footprints and spatial independence are unverified.')
  writeLines(caption, file.path(out, 'Study_area_caption.txt'))
  writeLines(c('# Study-area map revision', '', caption, '', '## Audit and remaining verification', 'All 28 ledger coordinates were checked against original worksheet degree-minute headers. Agreement is transcription verification, not field verification. Ngangao_06 remains assigned to Ngangao, as confirmed by the user. Its coordinate is unchanged; the source coordinate-versus-description discrepancy is retained in coordinate_source_audit.csv and is not a separate map category. Source workbooks and analytical data were not edited.', 'Panel (c) restores the independent analysis/Datasets/mapping/Study sites.csv reference points; they are explicitly distinguished from the worksheet coordinates in panel (d). Their original coordinates are retained in study_site_coordinate_sources.csv.', 'See the coordinate source audit, panel extents, scale-bar verification and environmental layer metadata in tables/.', 'Only the retained multipanel PNG and PDF are generated.'), file.path(out, 'README.md'))
  capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
  message('Revised map and audits saved in ', out)
  invisible(TRUE)
}
if (sys.nframe() == 0L) {
  script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))
  start <- if (length(script) && file.exists(script[1])) dirname(script[1]) else getwd()
  root <- normalizePath(start, winslash = '/', mustWork = TRUE)
  while (!file.exists(file.path(root, 'analysis/scripts/study_area.R'))) {
    parent <- dirname(root)
    if (identical(parent, root)) stop('Cannot locate the Taita project from ', start, call. = FALSE)
    root <- parent
  }
  run_study_area(root)
}
