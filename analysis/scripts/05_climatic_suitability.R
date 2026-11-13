# ============================================================
# 0. Project paths and output files
# ============================================================

# Shared output layout (embedded so each question remains self-contained).
rq_dirs <- function(root, question) {
  out <- file.path(root, 'analysis', 'outputs', paste('research question', question))
  for (p in c(out, file.path(out, c('figures', 'tables'))))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  out
}
rq_relative <- function(name) {
  name <- gsub('\\', '/', name, fixed = TRUE)
  is_figure <- grepl('[.](pdf|png|svg|jpg|jpeg|tif|tiff)$', name, ignore.case = TRUE)
  is_table <- grepl('[.](csv|tsv|xlsx|xls)$', name, ignore.case = TRUE) |
  grepl('^Table_.*[.]md$', basename(name))
  ifelse(grepl('^(figures|tables)/', name), name,
    ifelse(is_figure, paste0('figures/', name),
      ifelse(is_table, paste0('tables/', basename(name)), name)))
}
rq_file <- function(out, ...) {
  base <- sub('/figures(/.*)?$', '', gsub('\\', '/', out, fixed = TRUE))
  path <- file.path(base, rq_relative(file.path(...)))
  for (p in unique(dirname(path)))dir.create(p, recursive = TRUE, showWarnings = FALSE)
  path
}
rq_input <- function(root, ...) {
  file.path(root, 'analysis/Datasets/audits', basename(file.path(...)))
}

# ============================================================
# 1. Hs bool
# ============================================================

# Research question 5: Where do exploratory Maxent models predict suitable environmental conditions for each focal species, and how sensitive are those predictions to the calibration area and spatial thinning of occurrence records?
# Reviewed field and GBIF occurrences; crossed Maxent and separate legacy cached-map mode.

hs_bool <- function(x) !is.na(x) & tolower(trimws(as.character(x))) %in% c('true', '1', 'yes')

# ============================================================
# 2. Hs read
# ============================================================

hs_read <- function(p) {
  lines <- readLines(p, encoding = 'UTF-8', warn = FALSE)
  if (length(lines))lines[1] <- sub('^\ufeff', '', lines[1])
  read.csv(text = paste(lines, collapse = '\n'), stringsAsFactors = FALSE, check.names = FALSE, na.strings = c('', 'NA'), encoding = 'UTF-8')
}

# ============================================================
# 3. Hs coord ok
# ============================================================

hs_coord_ok <- function(lon, lat)is.finite(lon) & is.finite(lat) & abs(lon) <= 180 & abs(lat) <= 90 & !(lon == 0 & lat == 0)

# ============================================================
# 4. Hs reason
# ============================================================

hs_reason <- function(...) {
  x <- list(...);
  n <- length(x[[1]])
  vapply(seq_len(n), function(i)paste(vapply(x, function(z)if (is.na(z[i]))'' else z[i], character(1))[vapply(x, function(z)!is.na(z[i]) && nzchar(z[i]), logical(1))], collapse = '; '), character(1))
}

# ============================================================
# 5. Run habitat suitability
# ============================================================

run_habitat_suitability <- function(root) {
  root <- normalizePath(root);
  here <- file.path(root, 'analysis/scripts')
  .libPaths(c(file.path(here, '.r-library'), .libPaths()))
  if (!requireNamespace('jsonlite', quietly = TRUE))stop('jsonlite is required')
  input <- file.path(root, 'analysis/Datasets/habitat_suitability');
  out <- rq_dirs(root, 5)
  dir.create(input, recursive = TRUE, showWarnings = FALSE);
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  seed_csv <- function(x, n) {p <- file.path(input, paste0(n, '.csv'));
    if (!file.exists(p))write.csv(x, p, row.names = FALSE, na = 'NA');
    p}
  config_path <- file.path(input, 'config.json')
  if (!file.exists(config_path)) {
    spec <- list(calibration_mask = NULL, calibration_rationale = '', predictor_rationale = '',
      spatial_design_rationale = '', bias_raster = NULL, bias_rationale = '', prediction_mask = NULL,
      predictors = c('bio1', 'bio12', 'bio15'), max_coordinate_uncertainty_m = 1000,
      min_unique_cells = 20, min_train_presences = 10, min_test_presences = 3,
      spatial_block_m = 20000, exclusion_buffer_m = 5000, folds = 4, background_n = 4000,
      thinning_m = c(0, 1000), feature_classes = c('l', 'lq'), regularisation = c(1, 2, 4),
      min_outer_auc = .7, max_outer_omission = .25)
    jsonlite::write_json(list(seed = 5162026, environment = 'analysis/Datasets/legacy_habitat/environment/baseline.tif',
        species = list(Coffea_fadenii = spec, Afrocarpus_usambarensis = spec)), config_path, auto_unbox = TRUE, pretty = TRUE, null = 'null')
  }
  cfg <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
  # Shared source-row audit, not legacy climate-model inclusion flags.
  bp <- file.path(root, 'analysis/Datasets/audits')
  if (!all(file.exists(file.path(bp, c('record_reconciliation.csv', 'larger_stem_row_audit.csv', 'juvenile_row_audit.csv'))))) {
    cached <- file.path(root, 'analysis/Datasets/legacy_habitat/combined_source_occurrences.csv')
    if (!file.exists(cached))stop('Neither reconciled inputs nor cached occurrence audit is available')
    hs_restoration_maps(root, out, hs_read(cached))
    writeLines(c('# Restoration mapping run', '',
        'Generated exploratory restoration screening maps in restoration_maps/. Reconciled input files are unavailable; the existing occurrence audit and checksummed climate rasters were used. No Maxent refitting or new occurrence verification was performed.'), rq_file(out, 'MAP_RUN_STATUS.md'))
    message('Restoration maps complete; see ', rq_file(out, 'restoration_maps'))
    return(invisible(NULL))
  }
  q <- hs_read(file.path(bp, 'record_reconciliation.csv'))
  a <- hs_read(file.path(bp, 'larger_stem_row_audit.csv'));
  j <- hs_read(file.path(bp, 'juvenile_row_audit.csv'))
  species <- c('Coffea fadenii', 'Afrocarpus usambarensis');
  field <- list();
  k <- 0L
  for (sp in species)for (i in seq_len(nrow(q))) {
    aa <- a[a$plot_id == q$plot_id[i] & a$analysis_name == tolower(sp) & hs_bool(a$woody_eligible),]
    jj <- j[j$plot_id == q$plot_id[i] & j$analysis_name == tolower(sp) & hs_bool(j$woody_eligible),]
    header <- !is.na(q$focal_species_header[i]) && tolower(q$focal_species_header[i]) == tolower(sp)
    if (!nrow(aa) && !nrow(jj) && !header)next
    k <- k + 1L;
    field[[k]] <- data.frame(audit_id = paste('FIELD', gsub(' ', '_', sp), q$plot_id[i], sep = ':'), species = sp,
      source = 'Field worksheet', source_id = q$plot_id[i], occurrence_id = NA_character_, event_date = '2021-06',
      lon = q$lon[i], lat = q$lat[i], uncertainty_m = q$positional_uncertainty_m[i],
      original_name = q$focal_species_header[i], cultivated = hs_bool(q$cultivated[i]),
      listed_evidence = nrow(aa)>0 || nrow(jj)>0, primary_context = hs_bool(q$primary_eligible[i]),
      source_flags = q$exclusion_reason[i], source_file = q$source_file[i],
      evidence = paste0('Larger-stem rows=', nrow(aa), '; juvenile count=', sum(jj$count), '; header=', header))
  }
  field <- do.call(rbind, field)
  gbif_path <- file.path(root, 'analysis/Datasets/legacy_habitat/source_data/occurrence_review_queue.csv')
  if (!file.exists(gbif_path))stop('Cached GBIF occurrence queue missing: ', gbif_path)
  g <- hs_read(gbif_path);
  g <- g[g$source == 'GBIF' & g$species %in% species,]
  gbif <- data.frame(audit_id = paste('GBIF', g$species, g$record_id, sep = ':'), species = g$species, source = 'GBIF',
    source_id = as.character(g$record_id), occurrence_id = g$occurrence_id, event_date = g$date,
    lon = as.numeric(g$longitude), lat = as.numeric(g$latitude), uncertainty_m = as.numeric(g$uncertainty_m),
    original_name = g$scientific_name, cultivated = hs_bool(g$cultivated), listed_evidence = TRUE, primary_context = TRUE,
    source_flags = paste(g$review_flags, g$issues, sep = '; '), source_file = g$url,
    evidence = paste('basis=', g$basis, 'institution=', g$institution, 'catalog=', g$catalog))
  raw <- rbind(field, gbif);
  stopifnot(!anyDuplicated(raw$audit_id));
  put(raw, 'combined_source_occurrences')
  review_template <- data.frame(audit_id = raw$audit_id, accepted_species = raw$species,
    reviewed_lon = raw$lon, reviewed_lat = raw$lat, reviewed_uncertainty_m = raw$uncertainty_m,
    taxonomy_verified = FALSE, wild_status = 'unresolved', coordinate_verified = FALSE,
    context_resolved = FALSE, include_for_modelling = FALSE, duplicate_of = NA_character_,
    review_evidence = NA_character_, reviewer = NA_character_, review_date = NA_character_)
  review_path <- seed_csv(review_template, 'occurrence_review')
  rv <- hs_read(review_path);
  stopifnot(!anyDuplicated(rv$audit_id))
  required <- names(review_template);
  if (!all(required %in% names(rv)))stop('Review schema incomplete')
  if (any(!rv$audit_id %in% raw$audit_id))stop('Review has unknown occurrence IDs; reconcile the source first')
  new_rows <- review_template[!review_template$audit_id %in% rv$audit_id,]
  if (nrow(new_rows)) {
    rv <- rbind(rv, new_rows)
    write.csv(rv, review_path, row.names = FALSE, na = 'NA')
  }
  at <- match(raw$audit_id, rv$audit_id);
  r <- rv[at,]
  audited <- cbind(raw, r[, setdiff(names(r), 'audit_id'), drop = FALSE])
  audited$coordinate_numeric_valid <- hs_coord_ok(audited$reviewed_lon, audited$reviewed_lat)
  supplied <- function(x)!is.na(x) & nzchar(trimws(x))
  audited$review_ready <- hs_bool(r$taxonomy_verified) & r$wild_status == 'wild' & hs_bool(r$coordinate_verified) &
  hs_bool(r$context_resolved) & hs_bool(r$include_for_modelling) & supplied(r$review_evidence) & supplied(r$reviewer) & supplied(r$review_date)
  audited$review_ready[is.na(audited$review_ready)] <- FALSE
  audited$reason <- hs_reason(ifelse(hs_bool(r$taxonomy_verified), '', 'Taxonomic verification pending'),
    ifelse(!is.na(r$wild_status) & r$wild_status == 'wild', '', 'Wild context unresolved or not wild'),
    ifelse(hs_bool(r$coordinate_verified), '', 'Coordinates unverified'),
    ifelse(audited$coordinate_numeric_valid, '', 'Invalid coordinates'),
    ifelse(audited$review_ready, '', 'Review evidence, context or inclusion unresolved'),
    ifelse(raw$cultivated, 'Source records cultivation; must resolve explicitly', ''))
  # A cultivated source can only be corrected through explicit review, never silently.
  audited$review_ready <- audited$review_ready & !raw$cultivated
  audited$review_ready <- audited$review_ready & audited$coordinate_numeric_valid & !supplied(r$duplicate_of) &
  r$accepted_species == raw$species
  audited$review_ready[is.na(audited$review_ready)] <- FALSE
  audited$reason <- hs_reason(audited$reason, ifelse(supplied(r$duplicate_of), 'Duplicate identified in review', ''),
    ifelse(!is.na(r$accepted_species) & r$accepted_species == raw$species, '', 'Taxon changed; reconcile source before inclusion'))
  put(audited, 'occurrence_review_audit');
  put(review_template, 'current_occurrence_review_template')
  # Blank field forms are inputs, not invented measurements.
  sites <- data.frame(site_id = character(), species = character(), lon = numeric(), lat = numeric(),
    coordinate_verified = logical(), canopy_assessment = character(), dry_season_moisture = character(),
    substrate = character(), competing_vegetation = character(), threats_and_actions = character(),
    connectivity = character(), land_access = character(), management_commitment = character(),
    planting_provenance = character(), field_assessment_complete = logical(), notes = character())
  sites_path <- seed_csv(sites, 'candidate_sites')
  independent_path <- seed_csv(data.frame(species = character(), survey_id = character(), lon = numeric(), lat = numeric(),
      observed_presence = integer(), coordinate_verified = logical(), survey_comparable = logical(),
      survey_date = character(), evidence = character()), 'independent_validation')
  trial_path <- seed_csv(data.frame(site_id = character(), species = character(), block = character(),
      treatment = character(), plant_id = character(), provenance = character(), planting_date = character(),
      census_date = character(), status = character(), height_cm = numeric(), diameter_mm = numeric(),
      reproductive_status = character(), notes = character()), 'pilot_planting_monitoring')
  source(file.path(here, 'support/habitat_suitability_engine.Rinc'), local = TRUE)
  statuses <- list();
  models <- list()
  for (sp in species) {
    stem <- gsub(' ', '_', sp);
    s <- cfg$species[[stem]]
    if (is.null(s))stop('Missing species configuration: ', sp)
    v <- audited[audited$species == sp,]
    u <- is.finite(v$reviewed_uncertainty_m) & v$reviewed_uncertainty_m >= 0 &
    v$reviewed_uncertainty_m <= s$max_coordinate_uncertainty_m
    v$eligible_before_environment <- v$review_ready & u
    v$reason <- hs_reason(v$reason, ifelse(u, '', 'Uncertainty missing or above configured tolerance'))
    v$eligible_before_environment[is.na(v$eligible_before_environment)] <- FALSE
    put(v, paste0(stem, '_eligibility'))
    reasons <- character()
    if (sum(v$eligible_before_environment)<s$min_unique_cells)reasons <- c(reasons, 'Too few reviewed occurrences before cell deduplication')
    for (n in c('calibration_rationale', 'predictor_rationale', 'spatial_design_rationale'))
    if (is.null(s[[n]]) || !nzchar(s[[n]]))reasons <- c(reasons, paste('Missing', n))
    if (is.null(s$calibration_mask) || !file.exists(file.path(root, s$calibration_mask)))reasons <- c(reasons, 'Documented species calibration mask not supplied')
    for (p in c('terra', 'sf', 'maxnet'))if (!requireNamespace(p, quietly = TRUE))reasons <- c(reasons, paste('Package unavailable:', p))
    if (!file.exists(file.path(root, cfg$environment)))reasons <- c(reasons, 'Predictor raster missing')
    if (length(reasons)) {
      statuses[[sp]] <- data.frame(species = sp, source_records = nrow(v), reviewed_eligible = sum(v$eligible_before_environment),
        status = 'Audit complete; predictive modelling withheld', reason = paste(reasons, collapse = '; '),
        independent_validation_complete = FALSE, demonstrated_restoration_suitability = FALSE)
      next
    }
    result <- tryCatch(hs_model_species(v[v$eligible_before_environment,], s, cfg, root, out, sp,
        hs_read(sites_path), hs_read(independent_path)), error = function(e)list(status = paste('Model withheld:', conditionMessage(e)), independent = FALSE))
    statuses[[sp]] <- data.frame(species = sp, source_records = nrow(v), reviewed_eligible = sum(v$eligible_before_environment),
      status = result$status, reason = if (is.null(result$reason))''else result$reason,
      independent_validation_complete = isTRUE(result$independent), demonstrated_restoration_suitability = FALSE)
  }
  status <- do.call(rbind, statuses);
  put(status, 'species_model_status')
  # Always render available exploratory evidence separately from reviewed Maxent outputs.
  map_folder <- hs_restoration_maps(root, out, raw)
  # Profiles include source uncertainty instead of promoting pending points to training data.
  if (requireNamespace('terra', quietly = TRUE) && file.exists(file.path(root, cfg$environment))) {
    env <- terra::rast(file.path(root, cfg$environment))
    valid <- hs_coord_ok(raw$lon, raw$lat)
    vals <- matrix(NA_real_, nrow(raw), terra::nlyr(env));
    colnames(vals) <- names(env)
    if (any(valid)) {
      points <- terra::vect(raw[valid,], geom = c('lon', 'lat'), crs = 'EPSG:4326')
      points <- terra::project(points, terra::crs(env))
      vals[valid,] <- as.matrix(terra::extract(env, points, ID = FALSE))
    }
    put(cbind(raw[, c('audit_id', 'species', 'source', 'uncertainty_m')], review_ready = audited$review_ready, vals), 'provisional_environmental_profiles')
  }
  writeLines(c('# Habitat suitability: audit and gated modelling', '',
      'Separate species models combine reviewed wild field and cached GBIF records. Source coordinates are never multiplied by plant counts. Legacy model eligibility is not treated as verification. No GBIF download is performed; the cached queue and its provenance are retained.',
      'Review file and species configuration: analysis/Datasets/habitat_suitability. Existing review files are not overwritten. Reconcile new source IDs using current_occurrence_review_template.csv. Unknown coordinate uncertainty is not zero; supply documented accuracy appropriate to the raster resolution.',
      'Source cultivation excludes primary modelling. Taxon corrections and cultivated-status disputes require upstream source reconciliation. Cross-source duplicate events must be marked duplicate_of in the review; same raster-cell occurrences are additionally collapsed. Field records with unresolved context require explicit evidence before inclusion.',
      'Calibration masks are species-specific rasters aligned exactly to the predictor grid (positive values inside, NA or zero outside). Their biological rationale, predictor rationale and spatial-design rationale are required. Prediction masks use the same grid and define the restoration decision extent. No arbitrary buffer is substituted for a calibration region.',
      'Default minimum 20 occupied cells, four geographic folds, training/test minima, 20 km blocks and 5 km exclusion buffers are configurable operational safeguards, not universal published adequacy thresholds. Audit environmental coverage and justify these settings separately for each species.',
      'Maxent is fitted with maxnet. Linear and linear-quadratic classes and regularisation 1, 2, 4 are tuned by inner geographic folds within each outer training partition. All predictors are prespecified; no full-data response-based screening. A fixed linear regularisation-1 model is the benchmark. Fold failures and warnings are retained; incomplete comparisons do not qualify for mapping.',
      'Bias sensitivities compare within-cell deduplication and configured distance thinning, with area-weighted uniform background and, only if a documented sampling-effort raster exists, effort-weighted background. Missing effort data is reported, not invented. Geographic blocks are allocated to contiguous coordinate clusters; training presences/background are buffered from held-out presences. Fold memberships and sizes are written.',
      'Held-out AUC is presence-background discrimination, not true-absence specificity. Omission uses a training-presence tenth-percentile threshold. Report every fold. Relative cloglog scores are not occurrence, establishment or survival probabilities. No PA probability calibration is claimed.',
      'Mapping requires complete outer comparisons, configured discrimination/omission checks and improvement over the linear benchmark. These checks are screening criteria, not proof of ecological validity. Predictor-range extrapolation is masked from candidate scores and mapped separately; this diagnostic does not detect all multivariate novelty. Between-variant spread describes modelling sensitivity, not a confidence interval.',
      'Independent field surveys must be comparable, documented and outside every modelling occurrence cell. They evaluate observed occurrence; non-detection is not detection-corrected absence. Candidate scores and field-feasibility columns remain separate; no weighted restoration score is fabricated.',
      'Pilot planting forms do not represent collected observations. Replicated planting blocks, provenance, survival, growth and eventual reproduction are needed to demonstrate restoration performance. Existing climate projections are not reused as validated restoration maps.',
      '', 'References:',
      'Kramer-Schadt et al. (2013), Diversity and Distributions. https://doi.org/10.1111/ddi.12096',
      'Finn et al. (2024), Journal of Applied Ecology. https://doi.org/10.1111/1365-2664.14616',
      'Diez et al. (2014), Journal of Ecology. https://doi.org/10.1111/1365-2745.12215',
      'Hunter-Ayad et al. (2020), Journal of Applied Ecology. https://doi.org/10.1111/1365-2664.13629',
      'Implementation: https://cran.r-project.org/package=maxnet',
      '', 'Run: Rscript --vanilla "analysis/scripts/05_climatic_suitability.R"'), rq_file(out, 'METHODS.md'))

  inputs <- c(file.path(here, '05_climatic_suitability.R'), file.path(here, 'support/habitat_suitability_engine.Rinc'),
    config_path, review_path, sites_path, independent_path, trial_path, gbif_path,
    file.path(bp, c('record_reconciliation.csv', 'larger_stem_row_audit.csv', 'juvenile_row_audit.csv')),
    file.path(root, cfg$environment))
  for (s in cfg$species)for (n in c('calibration_mask', 'prediction_mask', 'bias_raster'))if (!is.null(s[[n]]))inputs <- c(inputs, file.path(root, s[[n]]))
  inputs <- unique(inputs[file.exists(inputs)])
  put(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))), 'input_and_code_hashes')
  message('Habitat-suitability workflow complete: ', out)
  invisible(status)
}

# ============================================================
# 6. Hs restoration maps
# ============================================================

# Standalone rendering of available exploratory climate evidence.
hs_restoration_maps <- function(root, out, occurrences) {
  if (!requireNamespace('terra', quietly = TRUE))stop('terra is required for maps')
  src <- file.path(root, 'analysis/Datasets/legacy_habitat')
  folder <- rq_file(out, 'restoration_maps');
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  manifest <- hs_read(file.path(src, 'raster_manifest.csv'))
  metrics <- hs_read(file.path(src, 'area_change_all_scenarios.csv'))
  q <- hs_read(rq_input(root, 'analysis/population_distribution/record_reconciliation.csv'))
  q <- q[hs_bool(q$primary_eligible) & hs_coord_ok(q$lon, q$lat),]
  window <- terra::ext(min(q$lon) - .08, max(q$lon) + .08, min(q$lat) - .08, max(q$lat) + .08)
  sources <- character();
  products <- character();
  models <- list();
  candidates <- list();
  summaries <- list()
  verify <- function(n) {
    path <- file.path(src, 'rasters', n);
    i <- match(n, manifest$file)
    if (is.na(i) || !file.exists(path) || unname(tools::md5sum(path)) != manifest$md5[i])stop('Missing or changed raster: ', n)
    sources <<- c(sources, path);
    path
  }
  for (sp in c('Coffea fadenii', 'Afrocarpus usambarensis')) {
    stem <- gsub(' ', '_', sp);
    full <- terra::rast(verify(paste0(stem, '_baseline_score.tif')))
    score <- terra::crop(full, window);
    sd <- terra::crop(terra::rast(verify(paste0(stem, '_small_model_sd.tif'))), score)
    cover <- terra::crop(terra::rast(verify('aligned_tree_cover_percent.tif')), score)
    terra::compareGeom(score, sd, cover, stopOnError = TRUE)
    th <- vapply(c('q05', 'q10'), function(t) {v <- unique(metrics$threshold_value[metrics$species == sp & metrics$threshold == t]);
        stopifnot(length(v) == 1, is.finite(v));
        v}, numeric(1))
    stopifnot(th[1] <= th[2]);
    v <- terra::values(score, mat = FALSE)
    stopifnot(all(v[is.finite(v)] >= 0 & v[is.finite(v)] <= 1))
    cl <- ifelse(is.finite(v), ifelse(v >= th[2], 3L, ifelse(v >= th[1], 2L, 1L)), NA_integer_)
    screen <- score;
    terra::values(screen) <- cl
    maps <- list(relative_climate_score = score, threshold_screen = screen, model_pair_sd = sd)
    for (n in names(maps)) {
      path <- rq_file(folder, paste0(stem, '_', n, '.tif'));
      terra::writeRaster(maps[[n]], path, overwrite = TRUE);
    products <- c(products, path)
    }
    ids <- which(cl == 3);
    xy <- terra::xyFromCell(score, ids);
    cv <- terra::values(cover, mat = FALSE);
    sv <- terra::values(sd, mat = FALSE)
    candidates[[sp]] <- data.frame(species = rep(sp, length(ids)), cell_id = terra::cellFromXY(full, xy),
      longitude = xy[, 1], latitude = xy[, 2], relative_climate_score = v[ids], model_pair_sd = sv[ids],
      circa2020_tree_cover_percent = cv[ids], threshold_q05 = rep(th[1], length(ids)), threshold_q10 = rep(th[2], length(ids)),
      field_assessed = rep(FALSE, length(ids)), establishment_tested = rep(FALSE, length(ids)))
    for (k in 1:3)summaries[[paste(sp, k)]] <- data.frame(species = sp, class = k, grid_cells = sum(cl == k, na.rm = TRUE), threshold_q05 = th[1], threshold_q10 = th[2])
    models[[sp]] <- list(full = full, score = score, screen = screen, sd = sd)
  }
  points_on_map <- function(sp, regional = FALSE) {
    z <- occurrences[occurrences$species == sp & hs_coord_ok(occurrences$lon, occurrences$lat),]
    if (!regional)z <- z[z$source == 'Field worksheet',]
    e <- if (regional)models[[sp]]$full else models[[sp]]$score
    z <- z[z$lon >= terra::xmin(e) & z$lon <= terra::xmax(e) & z$lat >= terra::ymin(e) & z$lat <= terra::ymax(e),]
    if (nrow(z))points(z$lon, z$lat, pch = ifelse(z$cultivated, 17, ifelse(z$source == 'GBIF', 1, 21)),
      bg = 'white', col = ifelse(z$cultivated, '#D55E00', 'black'), cex = if (regional).5 else .8)
  }
  foot <- function(title, note) {mtext(title, side = 3, outer = TRUE, line = 2, cex = 1.1, font = 2);
    mtext(note, side = 1, outer = TRUE, line = 1, cex = .7)}
  main <- function() {
    par(mfrow = c(2, 2), oma = c(7, 0, 4, 0), mar = c(4, 4, 3, 4))
    for (sp in names(models)) {
      terra::plot(models[[sp]]$score, col = hcl.colors(60, 'YlGnBu', rev = TRUE), range = c(0, 1),
        main = bquote(atop(italic(.(sp)), 'Relative climate score')), cex.main = .9, xlab = 'Longitude', ylab = 'Latitude');
      points_on_map(sp)
      terra::plot(models[[sp]]$screen, col = c('#E0E0E0', '#E9C46A', '#218C74'), breaks = c(.5, 1.5, 2.5, 3.5), legend = FALSE,
        main = bquote(atop(italic(.(sp)), 'Candidate-area threshold screen')), cex.main = .9, xlab = 'Longitude', ylab = 'Latitude');
      points_on_map(sp)
    }
    foot('Exploratory restoration candidate-area screening: Taita',
      'Screen: grey = below q05; gold = q05 only; green = above q05 and q10.\nCircles: provisional field records; orange triangles: cultivated. White raster cells: unresolved/no data.\nCached climate models; not validated planting sites. Extent: field-record envelope + 0.08 degrees; no forest boundaries inferred.')
  }
  sensitivity <- function() {
    par(mfrow = c(1, 2), oma = c(7, 0, 4, 0), mar = c(4, 4, 3, 4))
    for (sp in names(models)) {terra::plot(models[[sp]]$sd, col = hcl.colors(60, 'YlOrRd', rev = TRUE), main = bquote(italic(.(sp))), xlab = 'Longitude', ylab = 'Latitude');
    points_on_map(sp)}
    foot('Sensitivity among cached bivariate climate models',
      'SD across model pairs is not a confidence interval. Low SD does not establish accuracy.\nOccurrence verification, geographic transfer and planting performance remain unresolved.')
  }
  regional <- function() {
    par(mfrow = c(1, 2), oma = c(7, 0, 4, 0), mar = c(4, 4, 3, 4))
    for (sp in names(models)) {
      r <- terra::aggregate(models[[sp]]$full, fact = 4, fun = mean, na.rm = TRUE)
      terra::plot(r, col = hcl.colors(60, 'YlGnBu', rev = TRUE), range = c(0, 1), main = bquote(italic(.(sp))), xlab = 'Longitude', ylab = 'Latitude');
      points_on_map(sp, TRUE)
      rect(terra::xmin(window), terra::ymin(window), terra::xmax(window), terra::ymax(window), border = '#C51B7D', lwd = 2)
    }
    foot('Regional climate and occurrence context',
      'Hollow circles: GBIF; white-filled circles: field records; orange triangles: cultivated; pink box: Taita detail.\nPoints are not all verified or training records. Display raster aggregates 4 x 4 native cells.')
  }
  figures <- list(restoration_candidate_screening_maps = main, restoration_model_sensitivity_maps = sensitivity, regional_occurrence_context_maps = regional)
  for (n in names(figures))for (ext in c('pdf', 'png')) {
    path <- rq_file(folder, paste0(n, '.', ext));
    height <- if (n == 'restoration_candidate_screening_maps')10 else 7
    if (ext == 'pdf')pdf(path, width = 12, height = height)else png(path, width = 2400, height = height * 200, res = 200)
    tryCatch(figures[[n]](), finally = dev.off());
  products <- c(products, path)
  }
  write.csv(do.call(rbind, candidates), rq_file(folder, 'candidate_climate_grid_cells.csv'), row.names = FALSE)
  write.csv(do.call(rbind, summaries), rq_file(folder, 'screening_class_summary.csv'), row.names = FALSE)
  for (n in c('model_status.csv', 'geographic_validation.csv'))file.copy(file.path(src, n), rq_file(folder, paste0('source_', n)), overwrite = TRUE)
  sources <- unique(c(sources, file.path(src, c('raster_manifest.csv', 'area_change_all_scenarios.csv', 'model_status.csv', 'geographic_validation.csv')),
      file.path(root, 'analysis/scripts/05_climatic_suitability.R'), file.path(root, 'analysis/Datasets/legacy_habitat/combined_source_occurrences.csv'), rq_input(root, 'analysis/population_distribution/record_reconciliation.csv')))
  write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))), rq_file(folder, 'map_input_hashes.csv'), row.names = FALSE)
  write.csv(data.frame(path = products, md5 = unname(tools::md5sum(products))), rq_file(folder, 'map_output_manifest.csv'), row.names = FALSE)
  writeLines(c('# Restoration screening maps', '',
      'These newly rendered maps reuse existing ridge-penalised bivariate logistic climate models, NOT newly fitted Maxent models. Occurrence approvals remain unresolved and geographic validation has high omission and/or missing folds. See source_model_status.csv and source_geographic_validation.csv.',
      'Taita extent is the primary record envelope plus 0.08 degrees for display, not a calibration region. Native CHELSA 30-arc-second cells cannot identify planting microsites. Source masks omit water, missing inputs and baseline univariate extrapolation. White cells are not absence.',
      'q05 and q10 are fifth/tenth percentiles of historical training-presence scores, not establishment cutoffs. Classes show sensitivity to these two thresholds. candidate_climate_grid_cells.csv contains all cells above q10, not selected or field-verified planting sites.',
      'Tree cover is descriptive circa-2020 cover and is not a habitat exclusion mask. Model-pair SD is not a confidence interval or complete uncertainty. Regional display uses 4 x 4 aggregation; Taita outputs use the native grid.',
      'Circles are provisional records and can overlap; positional accuracy is unresolved. GBIF points are shown for context and are not necessarily model training points. Cultivated records are separately marked.',
      'Field feasibility assessments and replicated establishment trials are required before recommending planting sites.'), rq_file(folder, 'MAP_GUIDE.md'))
  message('Restoration screening maps generated: ', folder)
  invisible(folder)
}

# ============================================================
# 7. Mx auc
# ============================================================

# Exploratory Maxent helpers. All review flags remain unchanged.
mx_auc <- function(p, b) {
  if (!length(p) || !length(b) || any(!is.finite(c(p, b))))return(NA_real_)
  (sum(rank(c(p, b), ties.method = 'average')[seq_along(p)]) - length(p) * (length(p) + 1) / 2) / (length(p) * length(b))
}

# ============================================================
# 8. Mx distance
# ============================================================

mx_distance <- function(lon, lat, x, y) {
  r <- pi / 180
  2 * 6371000 * asin(pmin(1, sqrt(sin((lat - y) * r / 2)^2 + cos(lat * r) * cos(y * r) * sin((lon - x) * r / 2)^2)))
}

# ============================================================
# 9. Mx thin
# ============================================================

mx_thin <- function(p, d) {
  keep <- integer()
  for (i in seq_len(nrow(p)))if (!length(keep) || all(mx_distance(p$lon[keep], p$lat[keep], p$lon[i], p$lat[i]) >= d))keep <- c(keep, i)
  keep
}

# ============================================================
# 10. Mx far
# ============================================================

mx_far <- function(p, test, d) {
  if (!nrow(p))return(logical())
  result <- rep(TRUE, nrow(p))
  for (i in seq_len(nrow(test)))result <- result & mx_distance(p$lon, p$lat, test$lon[i], test$lat[i]) >= d
  result
}

# ============================================================
# 11. Mx fit
# ============================================================

mx_fit <- function(p, b, vars, classes, reg) {
  x <- rbind(p[, vars, drop = FALSE], b[, vars, drop = FALSE]);
  y <- c(rep(1, nrow(p)), rep(0, nrow(b)))
  if (any(vapply(x, function(v)length(unique(v))<2, logical(1))))stop('Constant training predictor')
  warnings <- character()
  m <- withCallingHandlers(maxnet::maxnet(y, x, maxnet::maxnet.formula(y, x, classes = classes), regmult = reg),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w));
      invokeRestart('muffleWarning')})
  if (!all(is.finite(m$betas)) || (!is.null(m$glmnet.fit$jerr) && m$glmnet.fit$jerr != 0))stop('Invalid fitted coefficients or glmnet error')
  list(model = m, vars = vars, lower = vapply(x, min, numeric(1)), upper = vapply(x, max, numeric(1)), warnings = paste(unique(warnings), collapse = '; '))
}

# ============================================================
# 12. Mx predict
# ============================================================

mx_predict <- function(m, x)as.numeric(predict(m$model, x[, m$vars, drop = FALSE], type = 'cloglog', clamp = FALSE))

# ============================================================
# 13. Mx novel
# ============================================================

mx_novel <- function(m, x) {
  xx <- as.matrix(x[, m$vars, drop = FALSE]);
  rowSums(sweep(xx, 2, m$lower, '<') | sweep(xx, 2, m$upper, '>'))>0
}

# ============================================================
# 14. Hs crossed inputs
# ============================================================

# Calibration masks and the prediction mask must share the predictor grid.
hs_crossed_inputs <- function(screened = FALSE) {
  stems <- c('Coffea_fadenii', 'Afrocarpus_usambarensis')
  occurrence <- if (screened)'screened_occurrences.csv'else 'verified_occurrences.csv'
  c(occurrence, 'crossed_design.csv', 'baseline.tif', 'prediction_mask.tif',
    unlist(lapply(stems, function(x)paste0(x, c('_narrow.tif', '_broad.tif')))))
}

# ============================================================
# 15. Run maxent exploratory
# ============================================================

run_maxent_exploratory <- function(root, input_dir = Sys.getenv('TAITA_HABITAT_INPUT_DIR',
    file.path(root, 'analysis/Datasets/habitat_suitability')), screened = FALSE) {
  out <- rq_dirs(root, 5)
  put <- function(x, name)write.csv(x, rq_file(out, paste0(name, '.csv')), row.names = FALSE, na = 'NA')
  species <- c('Coffea fadenii', 'Afrocarpus usambarensis')
  stems <- gsub(' ', '_', species)
  files <- hs_crossed_inputs(screened)
  occurrence_file <- files[1]
  optional <- c('forest_reference_points.csv', 'independent_assessment_sites.csv')
  paths <- setNames(file.path(input_dir, c(files, optional)), c(files, optional))
  missing <- files[!file.exists(paths[files])]
  if (length(missing)) stop("Required Question 5 inputs missing: ", paste(missing, collapse = ", "), call. = FALSE)
  .libPaths(c(file.path(root, 'analysis/scripts/.r-library'), .libPaths()))
  for (pkg in c('terra', 'maxnet'))if (!requireNamespace(pkg, quietly = TRUE))
  stop('Required package unavailable: ', pkg)
  require_columns <- function(x, columns) {
    if (!all(columns %in% names(x)))stop('Missing columns: ', paste(setdiff(columns, names(x)), collapse = ', '))
  }
  populated <- function(x)!is.na(x) & nzchar(trimws(x))
  design <- hs_read(paths[['crossed_design.csv']])
  require_columns(design, c('species', 'thinning_m', 'thinning_rationale', 'narrow_rationale',
      'broad_rationale', 'prediction_rationale'))
  stopifnot(!anyDuplicated(design$species), setequal(design$species, species),
    all(is.finite(design$thinning_m) & design$thinning_m>0),
    all(vapply(design[c('thinning_rationale', 'narrow_rationale', 'broad_rationale',
            'prediction_rationale')], function(x)all(populated(x)), logical(1))))
  if (!'stronger_thinning_m' %in% names(design))design$stronger_thinning_m <- NA_real_
  if (!'stronger_thinning_rationale' %in% names(design))design$stronger_thinning_rationale <- NA_character_
  stronger <- !is.na(design$stronger_thinning_m)
  stopifnot(all(is.finite(design$stronger_thinning_m[stronger])),
    all(design$stronger_thinning_m[stronger]>design$thinning_m[stronger]),
    all(populated(design$stronger_thinning_rationale[stronger])))
  vars <- c('bio1', 'bio12', 'bio15')
  env <- terra::rast(paths[['baseline.tif']])[[vars]]
  stopifnot(terra::is.lonlat(env), terra::same.crs(env, 'EPSG:4326'))
  template <- env[[1]]
  mask <- function(name) {
    x <- terra::rast(paths[[name]])
    stopifnot(terra::nlyr(x) == 1, terra::compareGeom(template, x, stopOnError = TRUE))
    which(terra::values(x, mat = FALSE) == 1)
  }
  values <- as.data.frame(terra::values(env))
  valid <- complete.cases(values) & Reduce('&', lapply(values, is.finite))
  grid <- mask('prediction_mask.tif')
  grid <- grid[valid[grid]]
  if (length(grid)<10)stop('Too few valid prediction cells')
  pred <- values[grid,, drop = FALSE]
  areas <- terra::values(terra::cellSize(template, unit = 'km'), mat = FALSE)
  cell_xy <- function(ids)if (length(ids))terra::xyFromCell(template, ids)else matrix(numeric(), 0, 2)
  xy_grid <- cell_xy(grid)
  half <- terra::res(template) / 2
  extent <- terra::ext(min(xy_grid[, 1]) - half[1], max(xy_grid[, 1]) + half[1],
    min(xy_grid[, 2]) - half[2], max(xy_grid[, 2]) + half[2])
  r <- hs_read(paths[[occurrence_file]])
  require_columns(r, c('audit_id', 'species', 'source', 'source_id', 'lon', 'lat', 'uncertainty_m',
      'taxonomy_verified', 'coordinate_verified', 'wild_verified', 'duplicate_of'))
  stopifnot(!anyDuplicated(r$audit_id), all(populated(r$audit_id)), all(populated(r$source_id)),
    all(r$species %in% species), all(r$source %in% c('Field', 'GBIF')))
  r$cell <- terra::cellFromXY(template, as.matrix(r[c('lon', 'lat')]))
  # Reserve assessment cells from both sources before deduplication.
  assessment <- paths[['independent_assessment_sites.csv']]
  assessment_cells <- integer()
  independent <- NULL
  if (file.exists(assessment)) {
    sites <- hs_read(assessment)
    require_columns(sites, c('label', 'lon', 'lat'))
    stopifnot(all(hs_coord_ok(sites$lon, sites$lat)))
    assessment_cells <- terra::cellFromXY(template, as.matrix(sites[c('lon', 'lat')]))
    if (all(c('species', 'verified_presence') %in% names(sites))) {
      check <- hs_bool(sites$verified_presence) & sites$species %in% species &
      !is.na(assessment_cells) & valid[assessment_cells]
      independent <- unique(data.frame(species = sites$species[check], cell = assessment_cells[check]))
    }
    assessment_cells <- assessment_cells[!is.na(assessment_cells)]
  }
  if (screened)require_columns(r, 'screen_pass')
  accepted <- if (screened)hs_bool(r$screen_pass) else
  hs_bool(r$taxonomy_verified) & hs_bool(r$coordinate_verified) & hs_bool(r$wild_verified)
  uncertainty_ok <- if (screened)is.na(r$uncertainty_m) else rep(FALSE, nrow(r))
  r$status <- 'Excluded: verification or screening, coordinates, uncertainty or duplicate flag'
  for (sp in species) {
    d <- design[design$species == sp,]
    ids <- which(r$species == sp & hs_coord_ok(r$lon, r$lat) &
      (uncertainty_ok | is.finite(r$uncertainty_m) & r$uncertainty_m >= 0 & r$uncertainty_m <= d$thinning_m / 2) &
      accepted &
      !populated(r$duplicate_of))
    r$status[ids] <- 'Excluded: missing predictors or outside raster'
    ids <- ids[which(valid[r$cell[ids]])]
    reserved <- r$cell[ids] %in% assessment_cells
    r$status[ids[reserved]] <- 'Excluded: independent assessment cell'
    ids <- ids[!reserved]
    ids <- ids[order(r$uncertainty_m[ids], r$source[ids] != 'Field', r$audit_id[ids])]
    source_keys <- paste(r$source[ids], r$source_id[ids], sep = '|')
    repeats <- duplicated(source_keys)
    r$status[ids[repeats]] <- 'Excluded: repeated source identifier'
    ids <- ids[!repeats]
    duplicate <- duplicated(r$cell[ids])
    r$status[ids] <- ifelse(duplicate, 'Excluded: duplicate predictor cell', 'Retained: baseline')
  }
  put(r, 'occurrence_inclusion_audit')
  put(as.data.frame(table(r$species, r$source, r$status),
      responseName = 'records'), 'occurrence_screening_counts')
  if (any(table(factor(r$species[r$status == 'Retained: baseline'], levels = species))<6))
  stop('Fewer than six eligible baseline cells for at least one species. See occurrence audit.')
  set.seed(5162026)
  models <- maps <- support <- validation <- sensitivity <- run_sensitivity <- candidates <- membership <- list()
  scenario_design <- cell_support <- scenario_area <- interactions <- independent_checks <- projection_summary <- list()
  percentile <- function(x)rank(x, ties.method = 'average') / length(x)
  grid_area <- areas[grid]
  metrics <- function(a, b, flag_a, flag_b) {
    union_area <- sum(grid_area[flag_a | flag_b])
    c(spearman = if (sd(a)>0 && sd(b)>0)cor(a, b, method = 'spearman')else NA_real_,
      jaccard = if (union_area>0)sum(grid_area[flag_a & flag_b]) / union_area else NA_real_,
      first_area_km2 = sum(grid_area[flag_a]), second_area_km2 = sum(grid_area[flag_b]),
      shared_km2 = sum(grid_area[flag_a & flag_b]), first_only_km2 = sum(grid_area[flag_a & !flag_b]),
      second_only_km2 = sum(grid_area[!flag_a & flag_b]))
  }
  thin <- function(x, d, seed) {
    set.seed(seed)
    y <- x[sample.int(nrow(x)),, drop = FALSE]
    y[mx_thin(y, d),, drop = FALSE]
  }
  for (i in seq_along(species)) {
    sp <- species[i];
    stem <- stems[i]
    spec <- design[design$species == sp,]
    distances <- c(baseline = 0, moderate = spec$thinning_m, stronger = spec$stronger_thinning_m)
    distances <- distances[is.finite(distances)]
    p <- r[r$species == sp & r$status == 'Retained: baseline',]
    eligible <- r[r$species == sp & r$status %in% c('Retained: baseline', 'Excluded: duplicate predictor cell'),]
    cell_support[[sp]] <- data.frame(species = sp, cell = p$cell, field = p$cell %in% eligible$cell[eligible$source == 'Field'],
      gbif = p$cell %in% eligible$cell[eligible$source == 'GBIF'])
    centres <- cell_xy(p$cell)
    p$lon <- centres[, 1];
    p$lat <- centres[, 2]
    # Test design support before fitting. A treatment is kept only if all repeats meet the fit minimum.
    seed_for <- function(run)61000L + 100L * i + run
    supported <- vapply(distances, function(distance)distance == 0 ||
      all(vapply(1:5, function(run)nrow(thin(p, distance, seed_for(run))) >= 6, logical(1))), logical(1))
    scenario_design[[sp]] <- data.frame(species = sp, treatment = names(distances),
      distance_m = unname(distances), supported = unname(supported))
    put(do.call(rbind, scenario_design), 'crossed_scenarios')
    distances <- distances[supported]
    scenarios <- expand.grid(treatment = names(distances), area = c('narrow', 'broad'), stringsAsFactors = FALSE)
    scenarios$name <- paste(scenarios$area, scenarios$treatment, sep = '_')
    scenarios$distance <- unname(distances[scenarios$treatment])
    scenarios$runs <- ifelse(scenarios$distance>0, 5L, 1L)
    # Keep the evaluation buffer identical across all treatments.
    d <- max(c(spec$thinning_m, distances))
    p <- cbind(p, values[p$cell,, drop = FALSE])
    domains <- setNames(lapply(c('narrow', 'broad'), function(a)
        mask(paste0(stem, '_', a, '.tif'))), c('narrow', 'broad'))
    if (!all(p$cell %in% domains$narrow) || !all(domains$narrow %in% domains$broad) ||
      length(domains$narrow) >= length(domains$broad))
    stop('Calibration masks must be nested, distinct, and both contain all baseline cells: ', sp)
    breaks <- unique(as.numeric(quantile(p$lat, c(1 / 3, 2 / 3))))
    fold_of <- function(y)if (length(breaks) == 2)findInterval(y, breaks) + 1L else rep(NA_integer_, length(y))
    p$fold <- fold_of(p$lat)
    bg <- lapply(domains, function(ids) {
        ids <- ids[valid[ids]]
        if (length(ids)<100)stop('Too few usable background cells: ', sp)
        ids <- ids[sample.int(length(ids), min(2500, length(ids)), prob = areas[ids])]
        xy <- cell_xy(ids)
        cbind(data.frame(lon = xy[, 1], lat = xy[, 2], fold = fold_of(xy[, 2])), values[ids,, drop = FALSE])
      })
    scores <- novelties <- flags_q10 <- flags_q05 <- list()
    for (j in seq_len(nrow(scenarios))) {
      scenario <- scenarios[j,]
      b <- bg[[scenario$area]]
      for (run in seq_len(scenario$runs)) {
        seed <- seed_for(run)
        train <- if (scenario$distance>0)thin(p, scenario$distance, seed)else p
        id <- paste(stem, scenario$name, run, sep = '_')
        fit <- if (nrow(train)<6) 'Fewer than six training cells' else
        tryCatch(mx_fit(train, b, vars, 'l', 2), error = function(e)conditionMessage(e))
        support[[id]] <- data.frame(species = sp, scenario = scenario$name, run, seed, thinning_m = scenario$distance,
          baseline_cells = nrow(p), training_cells = nrow(train), field_cells = sum(train$source == "Field"),
          gbif_cells = sum(train$source == "GBIF"), background_cells = nrow(b),
          features = 'linear', regularisation = 2, error = if (is.character(fit))fit else '',
          threshold_training_q05 = NA_real_, threshold_training_q10 = NA_real_, fit_warning = if (is.character(fit))''else fit$warnings)
        membership[[id]] <- data.frame(species = sp, scenario = scenario$name, run,
          audit_id = train$audit_id, cell = train$cell, source = train$source)
        if (!is.character(fit)) {
          models[[id]] <- fit
          score <- mx_predict(fit, pred)
          thresholds <- quantile(mx_predict(fit, train), c(.05, .10), names = FALSE)
          test_cells <- independent$cell[independent$species == sp]
          if (length(test_cells)) {
            check_score <- mx_predict(fit, values[test_cells,, drop = FALSE])
            independent_checks[[id]] <- data.frame(species = sp, scenario = scenario$name, run,
              heldout_cells = length(test_cells), omitted_q05 = sum(check_score<thresholds[1]),
              omitted_q10 = sum(check_score<thresholds[2]))
          }
          support[[id]]$threshold_training_q05 <- thresholds[1]
          support[[id]]$threshold_training_q10 <- thresholds[2]
          scores[[scenario$name]] <- cbind(scores[[scenario$name]], score)
          flags_q05[[scenario$name]] <- cbind(flags_q05[[scenario$name]], score >= thresholds[1])
          flags_q10[[scenario$name]] <- cbind(flags_q10[[scenario$name]], score >= thresholds[2])
          novelties[[scenario$name]] <- cbind(novelties[[scenario$name]], mx_novel(fit, pred))
        }
        for (f in 1:3) {
          test <- p[which(p$fold == f),]
          fold_train <- p[which(p$fold != f),]
          fold_train <- fold_train[mx_far(fold_train, test, d),]
          if (scenario$distance>0)fold_train <- thin(fold_train, scenario$distance, seed + f * 1000L)
          fold_bg <- b[which(b$fold != f),]
          fold_bg <- fold_bg[mx_far(fold_bg, test, d),]
          omitted <- threshold <- omitted_q05 <- threshold_q05 <- NA_real_
          error <- 'Insufficient spatially separated training, test or background cells'
          warning <- ''
          if (nrow(fold_train) >= 6 && nrow(test) >= 2 && nrow(fold_bg) >= 50) {
            error <- tryCatch({
                m <- mx_fit(fold_train, fold_bg, vars, 'l', 2)
                warning <- m$warnings
                thresholds <- quantile(mx_predict(m, fold_train), c(.05, .10), names = FALSE)
                threshold_q05 <- thresholds[1];
                threshold <- thresholds[2]
                test_score <- mx_predict(m, test)
                omitted <- sum(test_score<threshold)
                omitted_q05 <- sum(test_score<threshold_q05)
                ''
              }, error = function(e)conditionMessage(e))
          }
          validation[[paste(id, f)]] <- data.frame(species = sp, scenario = scenario$name, run, fold = f,
            train_cells = nrow(fold_train), test_cells = nrow(test), background_cells = nrow(fold_bg),
            omitted, omission_fraction = if (nrow(test))omitted / nrow(test)else NA_real_,
            threshold_training_q10 = threshold, omitted_q05, threshold_training_q05 = threshold_q05, warning, error,
            test_ids = paste(test$audit_id, collapse = '|'))
        }
      }
    }
    model_support <- do.call(rbind, support)
    folds <- do.call(rbind, validation)
    by_run <- split(folds, paste(folds$species, folds$scenario, folds$run, sep = '|'))
    key <- paste(model_support$species, model_support$scenario, model_support$run, sep = '|')
    model_support$successful_folds <- vapply(by_run[key], function(x)sum(is.finite(x$omitted)), integer(1))
    model_support$heldout_omission_fraction <- vapply(by_run[key], function(x)
      if (all(is.finite(x$omitted)))sum(x$omitted) / sum(x$test_cells)else NA_real_, numeric(1))
    put(model_support, 'model_support')
    spatial_validation <- data.frame(
      Species = folds$species,
      Scenario = gsub('_', ' ', folds$scenario, fixed = TRUE),
      Repeat = folds$run,
      Fold = folds$fold,
      `Training cells` = folds$train_cells,
      `Test cells` = folds$test_cells,
      `Background cells` = folds$background_cells,
      `Omitted Q10` = folds$omitted,
      `Omitted Q05` = folds$omitted_q05,
      `Fit status` = ifelse(is.finite(folds$omitted), 'Fitted',
        ifelse(startsWith(folds$error, 'Insufficient spatially separated'),
          'Insufficient separated data', 'Fit failed')),
      check.names = FALSE
     )
    validation_csv <- capture.output(write.csv(spatial_validation, row.names = FALSE,
      na = '__MISSING_OMISSION__'))
    writeLines(gsub('__MISSING_OMISSION__', '\u2014', validation_csv, fixed = TRUE, useBytes = TRUE),
      rq_file(out, 'Table_S7_Geographic_evaluation_by_withheld_latitude_band.csv'), useBytes = TRUE)
    put(do.call(rbind, membership), 'training_localities')
    if (!identical(names(scores), scenarios$name) ||
      !identical(unname(vapply(scores, ncol, integer(1))), scenarios$runs) ||
      any(!vapply(scores, function(x)all(is.finite(x)), logical(1))))
    stop('At least one required fit failed or produced invalid scores. See model support table.')
    ranks <- lapply(scores, function(x)apply(x, 2, percentile))
    combined <- vapply(ranks, rowMeans, numeric(length(grid)))
    if (any(apply(combined, 2, function(x)length(unique(x))<2)))
    stop('A scenario has constant predictions and cannot rank candidate cells: ', sp)
    high <- vapply(flags_q10, function(x)rowMeans(x)>.5, logical(length(grid)))
    high_q05 <- vapply(flags_q05, function(x)rowMeans(x)>.5, logical(length(grid)))
    agreement <- rowSums(high)
    novel <- Reduce('|', lapply(novelties, function(x)rowSums(x)>0))
    levels <- length(distances)
    comparisons <- cbind(seq_len(levels), seq_len(levels) + levels)
    if (levels>1) {
      pairs <- t(combn(seq_len(levels), 2))
      comparisons <- rbind(comparisons, pairs, pairs + levels)
    }
    differences <- list()
    for (k in seq_len(nrow(comparisons))) {
      a <- comparisons[k, 1];
      b <- comparisons[k, 2]
      label <- paste(stem, k, sep = '_')
      info <- data.frame(species = sp, first = scenarios$name[a], second = scenarios$name[b],
        comparison = if (k <= levels)'Calibration extent'else 'Thinning', common_cells = length(grid))
      result <- function(flags)metrics(combined[, a], combined[, b], flags[, a], flags[, b])
      sensitivity[[label]] <- rbind(cbind(info, threshold = 'q10', as.list(result(high))),
        cbind(info, threshold = 'q05', as.list(result(high_q05))))
      differences[[label]] <- 2L * high[, a] + high[, b]
      for (run in seq_len(max(ncol(ranks[[a]]), ncol(ranks[[b]])))) {
        ia <- min(run, ncol(ranks[[a]]));
      ib <- min(run, ncol(ranks[[b]]))
        run_sensitivity[[paste(label, run)]] <- cbind(info, run, as.list(metrics(
              ranks[[a]][, ia], ranks[[b]][, ib], flags_q10[[a]][, ia], flags_q10[[b]][, ib])))
      }
    }
    if (levels>1) {
      x <- do.call(rbind, sensitivity)
      x <- x[x$species == sp & x$comparison == 'Thinning',]
      x$treatment_pair <- paste(sub('^[^_]+_', '', x$first), sub('^[^_]+_', '', x$second), sep = ' vs ')
      interactions[[sp]] <- merge(x[startsWith(x$first, 'narrow'), c('threshold', 'treatment_pair', 'spearman', 'jaccard')],
        x[startsWith(x$first, 'broad'), c('threshold', 'treatment_pair', 'spearman', 'jaccard')],
        by = c('threshold', 'treatment_pair'), suffixes = c('_narrow', '_broad'))
      interactions[[sp]] <- transform(interactions[[sp]], species = sp,
        jaccard_narrow_minus_broad = jaccard_narrow - jaccard_broad,
        spearman_narrow_minus_broad = spearman_narrow - spearman_broad)
    }
    scenario_area[[sp]] <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(j)
        data.frame(species = sp, scenario = scenarios$name[j], run = seq_len(scenarios$runs[j]),
          q10_area_km2 = colSums(flags_q10[[j]] * grid_area), q05_area_km2 = colSums(flags_q05[[j]] * grid_area),
          q10_scenario_km2 = sum(grid_area[high[, j]]), q05_scenario_km2 = sum(grid_area[high_q05[, j]]))))
    make_raster <- function(x, name) {
      z <- template
      v <- rep(NA_real_, terra::ncell(z));
      v[grid] <- x
      terra::values(z) <- v
      z <- terra::crop(z, extent)
      terra::writeRaster(z, rq_file(out, paste0(stem, '_', name, '.tif')), overwrite = TRUE)
      z
    }
    reference <- scores$narrow_baseline[, 1]
    contrast_rows <- do.call(rbind, sensitivity[paste(stem, seq_len(nrow(comparisons)), sep = '_')])
    contrast_rows <- contrast_rows[contrast_rows$threshold == 'q10',]
    most <- which.min(replace(contrast_rows$jaccard, is.na(contrast_rows$jaccard), Inf))
    maps[[sp]] <- list(n_scenarios = nrow(scenarios), design = scenarios,
      difference = make_raster(differences[[most]], 'strongest_scenario_disagreement'),
      difference_names = unlist(contrast_rows[most, c('first', 'second')], use.names = FALSE),
      both_source = cell_xy(cell_support[[sp]]$cell[cell_support[[sp]]$field & cell_support[[sp]]$gbif]),
      reference = make_raster(reference, 'reference_suitability'),
      agreement = make_raster(agreement, 'scenario_agreement'),
      scenario = lapply(setNames(seq_len(nrow(scenarios)), scenarios$name), function(j)
        make_raster(combined[, j], paste0(scenarios$name[j], '_percentile_rank'))),
      reference_novel = cell_xy(grid[novelties$narrow_baseline[, 1]]),
      novel = cell_xy(grid[novel]), domains = domains, points = r[r$species == sp & r$status == 'Retained: baseline',])
    for (j in seq_len(nrow(scenarios))) {
      make_raster(as.integer(rowSums(novelties[[j]])>0),
        paste0(scenarios$name[j], '_outside_training_range'))
    }
    make_raster(as.integer(novel), 'outside_any_training_range')
    for (scenario_name in c('narrow_baseline', 'broad_baseline'))
      make_raster(as.integer(high[, scenario_name]), paste0(scenario_name, '_q10_suitable'))
    make_raster(as.integer(agreement == nrow(scenarios) & !novel), 'shortlist_q10_suitable')
    take <- which(agreement == nrow(scenarios) & !novel)
    projection_summary[[sp]] <- data.frame(species = sp, scenarios = nrow(scenarios),
      prediction_cells = length(grid), prediction_km2 = sum(grid_area),
      extrapolated_cells = sum(novel), extrapolated_km2 = sum(grid_area[novel]),
      all_agree_cells = sum(agreement == nrow(scenarios)),
      all_agree_km2 = sum(grid_area[agreement == nrow(scenarios)]),
      candidate_cells = length(take), candidate_km2 = sum(grid_area[take]))
    xy <- cell_xy(grid[take])
    candidates[[sp]] <- data.frame(species = rep(sp, length(take)), cell = grid[take],
      longitude = xy[, 1], latitude = xy[, 2], reference_suitability = reference[take],
      agreeing_scenarios = agreement[take])
  }
  table <- do.call(rbind, sensitivity)
  hs_write_rank_sensitivity(table, rq_file(out, 'Table_2_Sensitivity_of_Maxent_predictions.csv'))
  contrasts_by_choice <- split(table, list(table$species, table$threshold, table$comparison), drop = TRUE)
  most_sensitive <- lapply(contrasts_by_choice, function(x) {
      if (!any(is.finite(x$jaccard)))return(NULL)
      x[which.min(x$jaccard),, drop = FALSE]
    })
  if (any(lengths(most_sensitive)>0L))
  put(do.call(rbind, most_sensitive), 'lowest_overlap_by_modelling_choice')
  put(do.call(rbind, cell_support), 'occurrence_cell_source_support')
  training_table <- do.call(rbind, membership)
  retained <- do.call(rbind, support)
  put(data.frame(
    Species = retained$species,
    Scenario = gsub('_', ' ', retained$scenario, fixed = TRUE),
    Repeat = retained$run,
    `Thinning (km)` = retained$thinning_m / 1000,
    `Training cells` = retained$training_cells,
    `Field cells` = retained$field_cells,
    `GBIF cells` = retained$gbif_cells,
    check.names = FALSE
  ), 'Table_S6_Occupied_predictor_cells_retained_for_Maxent_fitting')
  checks <- lapply(split(training_table, list(training_table$species, training_table$run), drop = TRUE), function(x) {
      a <- x$cell[x$scenario == 'broad_moderate'];
      b <- x$cell[x$scenario == 'broad_stronger']
      if (!length(a) || !length(b))return(NULL)
      data.frame(species = x$species[1], run = x$run[1], moderate_cells = length(a),
        stronger_cells = length(b), shared_cells = length(intersect(a, b)), same_cell_set = setequal(a, b))
    })
  cell_check <- do.call(rbind, checks)
  if (!is.null(cell_check))put(cell_check, 'moderate_vs_stronger_training_cells')
  area_table <- do.call(rbind, scenario_area)
  repeat_table <- do.call(rbind, run_sensitivity)
  put(area_table, 'scenario_suitable_area_km2')
  if (length(independent_checks))put(do.call(rbind, independent_checks), 'independent_presence_omission')
  if (length(interactions))put(do.call(rbind, interactions), 'calibration_by_thinning_interaction')
  put(do.call(rbind, run_sensitivity), 'thinning_repeat_sensitivity')
  put(do.call(rbind, candidates), 'candidate_screening_grid_cells')
  put(do.call(rbind, projection_summary), 'prediction_coverage_and_agreement')

  hs_crossed_figures(out, maps, template, grid, paths, repeat_table, area_table, screened)
  existing <- paths[file.exists(paths)]
  put(data.frame(path = existing, md5 = unname(tools::md5sum(existing))), 'crossed_input_hashes')

  invisible(list(models = models, sensitivity = table, output = out))
}

# ============================================================
# 16. Hs main results
# ============================================================


# ============================================================
# 17. Hs crossed figures
# ============================================================

hs_write_rank_sensitivity <- function(table, path) {
  baseline <- table$first == 'narrow_baseline' & table$second == 'broad_baseline'
  moderate <- table$species == 'Afrocarpus usambarensis' & (
    (table$first == 'narrow_baseline' & table$second == 'narrow_moderate') |
    (table$first == 'broad_baseline' & table$second == 'broad_moderate'))
  d <- table[which(table$threshold == 'q10' & (baseline | moderate)), ]
  d <- d[order(match(d$species, c('Coffea fadenii', 'Afrocarpus usambarensis')),
    match(d$second, c('broad_baseline', 'narrow_moderate', 'broad_moderate'))), ]
  area <- function(x) formatC(x, format = 'f', digits = 0, big.mark = ',')
  result <- data.frame(
    Species = sub('^([[:alpha:]])[^ ]+ ', '\\1. ', d$species),
    'Scenario contrast' = paste(gsub('_', ' ', d$first), 'to', gsub('_', ' ', d$second)),
    Jaccard = sprintf('%.2f', d$jaccard),
    Spearman = sprintf('%.4f', d$spearman),
    'Suitable area' = paste(area(d$first_area_km2), 'to', area(d$second_area_km2)),
    check.names = FALSE)
  result[['Scenario contrast']] <- sub('^n', 'N', sub('^b', 'B', result[['Scenario contrast']]))
  lines <- capture.output(write.csv(result, row.names = FALSE))
  lines[1] <- '"Species","Scenario contrast","Jaccard","Spearman \u03c1","Suitable area (km\u00b2)"'
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}

hs_crossed_figures <- function(out, maps, template, grid, paths, repeat_table, area_table, screened = FALSE) {
  frame_colour <- 'grey50'
  # Low-to-high colours sampled at 10% intervals from the reference colour bar.
  suit <- grDevices::colorRampPalette(c(
    '#FCCCF5', '#FEB9BF', '#FAA889', '#E0954E', '#B28D2B', '#82822E',
    '#577543', '#2D695D', '#125462', '#073960', '#001357'
  ))(50)
  vote_palette <- colorRampPalette(c('#E6E6E6', '#E8AE54', '#74B7A5', '#006B59'))
  half <- terra::res(template) / 2
  hatch <- function(xy, density = 17, col = '#494949') {
    if (nrow(xy))rect(xy[, 1] - half[1], xy[, 2] - half[2], xy[, 1] + half[1], xy[, 2] + half[2],
      density = density, angle = 45, col = col, border = NA)
  }
  polygon <- function(cells) {
    z <- template
    v <- rep(NA_real_, terra::ncell(z));
    v[cells] <- 1
    terra::values(z) <- v
    terra::as.polygons(terra::trim(z), aggregate = TRUE)
  }
  sites <- function(name) {
    if (!file.exists(paths[[name]]))return(NULL)
    z <- hs_read(paths[[name]])
    stopifnot(all(c('label', 'lon', 'lat') %in% names(z)), all(hs_coord_ok(z$lon, z$lat)))
    z
  }
  refs <- sites('forest_reference_points.csv')
  independent <- sites('independent_assessment_sites.csv')
  titles <- list(expression(italic(Coffea ~ fadenii)), expression(italic(Afrocarpus ~ usambarensis)))
  scenario_label <- function(x) {
    x <- gsub('_', ' ', x)
    paste0(toupper(substr(x, 1, 1)), substring(x, 2))
  }
  draw_occurrences <- function() {
    coordinate_title_cex <- 15/12
    layout(rbind(c(1, 0, 2), c(3, 0, 4), c(5, 5, 5)),
      widths = c(1, 0, 1), heights = c(4, 4, 1))
    par(mar = c(2.55, 3.3, 3.2, .8), oma = c(0, 1.7, 0, 0),
      cex = 1, cex.axis = 14/12, mgp = c(2, .55, 0))
    prediction <- polygon(grid)
    strip <- function(i, panel, subtitle) {
      bounds <- par('usr')
      y <- grconvertY(grconvertY(bounds[4], 'user', 'inches') +
        c(.5/25.4, .135, .365, .49), 'inches', 'user')
      rect(bounds[1], y[1], bounds[2], y[4], col = '#E5E5E5', border = 'grey50', xpd = NA)
      text(mean(bounds[1:2]), mean(y[c(1, 4)]),
        bquote(bold(.(paste0('(', letters[panel], ')'))) ~ .(titles[[i]][[1]]) ~ .(subtitle)),
        cex = 14/12, xpd = NA)
    }
    points_on_map <- function(m, size = 1) {
      g <- m$points$source == 'GBIF'
      points(m$points$lon, m$points$lat, pch = ifelse(g, 1, 17),
        col = ifelse(g, '#333333', '#9B2C6C'), cex = size)
      if (nrow(m$both_source))points(m$both_source, pch = 23, bg = '#5AAE61', cex = 1.05 * size)
      if (!is.null(independent))points(independent$lon, independent$lat, pch = 8, col = '#147D64', cex = size)
    }
    for (i in seq_along(maps)) {
      par(mar = c(2.55, if (i == 1L) 3.3 else 3.05, 3.2, if (i == 1L) .55 else .8))
      m <- maps[[i]]
      broad <- polygon(m$domains$broad)
      bounds <- as.vector(terra::union(terra::ext(broad), terra::ext(prediction)))
      plot(NA, xlim = bounds[1:2], ylim = bounds[3:4], asp = 1, xlab = '', ylab = '', bty = 'n')
      box(col = frame_colour)
      if (i == 1L) latitude_x <- grconvertX(
        grconvertX(par('usr')[1], 'user', 'inches') - .34 - 4/25.4, 'inches', 'ndc')
      strip(i, i, '(calibration areas)')
      terra::plot(broad, add = TRUE, col = '#F4F4F4', border = '#555555')
      terra::plot(polygon(m$domains$narrow), add = TRUE, col = NA, border = '#0072B2', lwd = 1.6)
      terra::plot(prediction, add = TRUE, col = NA, border = '#D55E00', lwd = 2)
      points_on_map(m)
    }
    e <- if (screened)c(38.12, 38.65, - 3.62, - 3.12)else as.vector(terra::ext(prediction))
    detail <- terra::crop(prediction, terra::ext(e))
    for (i in seq_along(maps)) {
      par(mar = c(2.55, if (i == 1L) 3.3 else 3.05, 3.2, if (i == 1L) .55 else .8))
      par(mai = par('mai') + c(1, 0, -1, 0)/25.4)
      m <- maps[[i]]
      plot(NA, xlim = e[1:2], ylim = e[3:4], asp = 1, xlab = '', ylab = '', bty = 'n')
      box(col = frame_colour)
      strip(i, i + 2L, '(Taita Hills detail)')
      terra::plot(detail, add = TRUE, col = '#FAF5F1', border = '#D55E00')
      points_on_map(m, size = 1.4)
      if (!is.null(refs)) {
        points(refs$lon, refs$lat, pch = 4, cex = 1.5, lwd = 1.4)
        west <- refs$lon <= mean(refs$lon)
        label_x <- refs$lon + ifelse(west, -.028, .028)
        label_y <- refs$lat - .065
        segments(refs$lon, refs$lat - .012, label_x, label_y, col = 'grey45')
        text(label_x[west], label_y[west], refs$label[west], pos = 2, offset = .18, cex = 1)
        text(label_x[!west], label_y[!west], refs$label[!west], pos = 4, offset = .18, cex = 1)
      }
    }
    par(mar = rep(0, 4))
    plot.new()
    labels <- c('Field training', 'GBIF training', 'Both sources in cell',
      if (!is.null(independent))'Independent assessment', if (!is.null(refs))'Forest reference')
    text(.5, .89, expression(bold(Longitude~(degree))), cex = coordinate_title_cex)
    occurrence_key <- legend(.5, .52, labels, pch = c(17, 1, 23, if (!is.null(independent))8, if (!is.null(refs))4),
      pt.bg = c(NA, NA, '#5AAE61', rep(NA, length(labels) - 3)),
      col = c('#9B2C6C', '#333333', 'black', if (!is.null(independent))'#147D64', if (!is.null(refs))'black'),
      horiz = TRUE, xjust = .5, yjust = .5, bty = 'n', cex = 14/12, pt.cex = 2.6, x.intersp = 1.5)
    segments(.05, .35, .95, .35, col = 'grey65')
    legend(occurrence_key$rect$left, .18, c('Broad calibration', 'Narrow calibration', 'Common prediction area'),
      col = c('#555555', '#0072B2', '#D55E00'), lty = 1, lwd = 7, seg.len = 3,
      horiz = TRUE, xjust = 0, yjust = .5, bty = 'n', cex = 14/12,
      text.width = NA, x.intersp = .55)
    par(fig = c(0, 1, 0, 1), new = TRUE, oma = rep(0, 4))
    plot.new()
    text(grconvertX(latitude_x, 'ndc', 'user'), grconvertY(.56, 'ndc', 'user'),
      expression(bold(Latitude~(degree))), srt = 90, cex = coordinate_title_cex)
  }
  draw_raster <- function(raster, colours, breaks) {
    e <- as.vector(terra::ext(raster))
    index <- cut(terra::values(raster, mat = FALSE), breaks,
      labels = FALSE, include.lowest = TRUE)
    pixels <- colours[index]
    pixels[is.na(pixels)] <- '#FFFFFF00'
    rasterImage(as.raster(matrix(pixels, nrow = nrow(raster), byrow = TRUE)),
      e[1], e[3], e[2], e[4], interpolate = FALSE)
  }
  map_axes <- function(column, top_row) {
    par(xpd = NA)
    box(col = frame_colour)
    longitude_ticks <- seq(34, 42, 2)
    axis(1, at = longitude_ticks, labels = FALSE, col = frame_colour)
    tick_top <- grconvertY(par('usr')[3], 'user', 'inches') - .19
    text(longitude_ticks, grconvertY(tick_top, 'inches', 'user'),
      labels = longitude_ticks, adj = c(.5, 1), cex = 15/12)
    if (column == 1L) axis(2, cex.axis = 15/12, col = frame_colour)
  }

  map_heading <- function(letter, species, subtitle, left = FALSE, species_below = FALSE) {
    bounds <- par('usr')
    centre <- if (left) bounds[1] else mean(bounds[1:2])
    alignment <- if (left) 0 else .5
    frame_top <- grconvertY(bounds[4], 'user', 'inches')
    subtitle_y <- frame_top + .10
    strip_y <- grconvertY(frame_top + c(.5/25.4, .64), 'inches', 'user')
    rect(bounds[1], strip_y[1], bounds[2], strip_y[2], col = '#E5E5E5',
      border = frame_colour, lwd = 1, xpd = NA)
    heading <- bquote(bold(.(paste0('(', letter, ')'))) ~ .(species))
    if (species_below) {
      heading <- bquote(bold(.(paste0('(', letter, ')'))) ~ .(subtitle))
      subtitle <- bquote(group("(", .(species), ")"))
    }
    text(centre, grconvertY(subtitle_y, 'inches', 'user'), subtitle,
      adj = c(alignment, 0), cex = if (species_below) 15/12 else 16.6/12, font = 1, xpd = NA)
    text(centre, grconvertY(subtitle_y + .25, 'inches', 'user'),
      heading,
      adj = c(alignment, 0), cex = 16.6/12, xpd = NA)
  }
  draw_main <- function() {
    par(mar = rep(0, 4), oma = rep(0, 4), cex.axis = 15/12, xpd = NA)
    shortlist_colour <- '#6A235D'
    agreement_fractions <- sort(unique(unlist(lapply(maps, function(m)
      (0:m$n_scenarios) / m$n_scenarios))))
    agreement_breaks <- c(-.01, head(agreement_fractions, -1) + diff(agreement_fractions)/2, 1.01)
    agreement_colours <- vote_palette(101)[1 + round(100 * agreement_fractions)]
    agreement_labels <- paste0(round(100 * agreement_fractions, 1), '%')
    map_window <- function(raster, column, top_row) {
      e <- as.vector(terra::ext(raster))
      vertical <- if (top_row) c(.63, .94) else c(.18, .49)
      device <- par('din')
      vertical <- (vertical * 15.18 - if (top_row) 1.18 else .68) / device[2]
      width <- 1.15 * diff(vertical) * device[2] * diff(e[1:2])/diff(e[3:4])/device[1]
      vertical <- mean(vertical) + 1.05 * (vertical - mean(vertical))
      if (!top_row) {
        upper_bottom <- (15.18 * (.63 - .31 * .025) - 1.18) / device[2]
        tick_height <- max(strheight(as.character(seq(34, 42, 2)), units = 'inches', cex = 15/12))
        title_height <- strheight(expression(bold(Longitude~(degree))), units = 'inches', cex = 17/12)
        lower_top <- upper_bottom - (.19 + tick_height + 3/25.4 + title_height + 4/25.4 + .64)/device[2]
        vertical <- vertical + lower_top - vertical[2]
      }
      half_gap <- 1/25.4/device[1]
      horizontal <- if (column == 1L) .5 - half_gap - c(width, 0) else
        .5 + half_gap + c(0, width)
      par(fig = c(0, 1, 0, 1), plt = c(horizontal, vertical),
        new = !(column == 1L && top_row), xpd = NA)
      plot.new()
      padding <- .07 * diff(e[3:4]) / (diff(vertical)*device[2] - .14)
      x_limits <- range(e[1:2] + c(-padding, padding), 42 + padding)
      plot.window(xlim = x_limits,
        ylim = e[3:4] + c(-padding, padding), xaxs = 'i', yaxs = 'i')
    }
    for (i in seq_along(maps)) {
      m <- maps[[i]]
      map_window(m$reference, i, TRUE)
      if (i == 2L) baseline_frame <- par('plt')
      draw_raster(m$reference, suit, seq(0, 1, length.out = length(suit) + 1))
      map_axes(i, TRUE)
      map_heading(letters[i], titles[[i]][[1]], '(baseline suitability)')
      hatch(m$reference_novel, density = 8, col = '#777777')
      if (!is.null(refs)) {
        points(refs$lon, refs$lat, pch = 21, bg = 'white', cex = 1.15 + 2/9)
        text(mean(refs$lon), mean(refs$lat), 'Taita Hills', pos = 3, cex = 14/12)
      }
      map_window(m$agreement, i, FALSE)
      if (i == 2L) agreement_frame <- par('plt')
      draw_raster(m$agreement / m$n_scenarios, agreement_colours, agreement_breaks)
      map_axes(i, FALSE)
      agreement_heading <- sprintf('(%d scenarios)', m$n_scenarios)
      map_heading(letters[i+2L], titles[[i]][[1]], agreement_heading)
      hatch(m$novel, density = 8, col = '#777777')
      candidate <- which(terra::values(m$agreement, mat = FALSE) == m$n_scenarios)
      if (nrow(m$novel)) candidate <- setdiff(candidate, terra::cellFromXY(m$agreement, m$novel))
      if (length(candidate)) {
        mask <- m$agreement
        v <- rep(NA_real_, terra::ncell(mask))
        v[candidate] <- 1
        terra::values(mask) <- v
        outline <- terra::as.polygons(mask, aggregate = TRUE)
        terra::plot(outline, add = TRUE, col = NA, border = '#FAF5ED', lwd = 2.7)
        terra::plot(outline, add = TRUE, col = NA, border = shortlist_colour, lwd = 1.6)
      }
      if (!is.null(refs)) {
        points(refs$lon, refs$lat, pch = 21, bg = 'white', cex = .8 + 2/9)
        text(mean(refs$lon), mean(refs$lat), 'Taita Hills', pos = 3, cex = 10.8/12)
      }
    }
    gap <- 2/25.4/par('din')[1]
    legend_left <- baseline_frame[2] + gap
    par(fig = c(legend_left, legend_left + .10, baseline_frame[3:4]), new = TRUE)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
    scale_y <- seq(0, 1, length.out = length(suit) + 1)
    rect(0, head(scale_y, -1), .09, tail(scale_y, -1), col = suit, border = NA)
    rect(0, 0, .09, 1, border = 'black')
    text(.17, seq(0, 1, length.out = 6), sprintf('%.1f', seq(0, 1, .2)),
      adj = c(0, .5), cex = 15/12)
    legend_left <- agreement_frame[2] + gap
    par(fig = c(legend_left, min(1, legend_left + .2), agreement_frame[3:4]), new = TRUE)
    plot.new()
    legend('left', agreement_labels, pch = 22, pt.bg = agreement_colours,
      pt.cex = 3, bty = 'n', cex = 15/12, y.intersp = 1.2)
    par(fig = c(0, 1, 0, .09 * 15.18/par('din')[2]), new = TRUE)
    plot.new()
    bottom_labels <- c('Unanimously suitable and within training ranges in every fit',
      'Outside training ranges')
    key_left <- .5 - (.06 + max(strwidth(bottom_labels, cex = 16/12)))/2
    segments(key_left, .66, key_left + .045, .66, col = shortlist_colour, lwd = 8, lend = 'butt')
    rect(key_left, .10, key_left + .045, .38, density = 8, angle = 45,
      col = '#777777', border = 'black')
    text(key_left + .06, c(.66, .24), bottom_labels, adj = c(0, .5), cex = 16/12)
    par(fig = c(0, 1, 0, 1), new = TRUE)
    plot.new()
    frame_left <- 1 - baseline_frame[2]
    latitude_x <- frame_left - (.36 + 2/25.4 + strheight('Latitude', units = 'inches',
      cex = 17/12, font = 2)/2)/par('din')[1]
    latitude_y <- c(mean(baseline_frame[3:4]), mean(agreement_frame[3:4]))
    text(grconvertX(latitude_x, 'ndc', 'user'), grconvertY(latitude_y, 'ndc', 'user'),
      'Latitude (\u00b0)\n', srt = 90, font = 2, cex = 17/12)
    tick_height <- max(strheight(as.character(seq(34, 42, 2)), units = 'inches', cex = 15/12))
    longitude_y <- c(baseline_frame[3], agreement_frame[3]) -
      (.19 + tick_height + 3/25.4)/par('din')[2]
    text(grconvertX(.5, 'ndc', 'user'), grconvertY(longitude_y, 'ndc', 'user'),
      expression(bold(Longitude~(degree))), adj = c(.5, 1), cex = 17/12)
  }
  draw_scenarios <- function() {
    par(family = 'Arial', fg = 'black', col = 'black', col.axis = 'black', col.lab = 'black')
    par(omi = c(.16, .48, 0, .72), bg = 'white')
    coffea <- maps[['Coffea fadenii']]$scenario
    afrocarpus <- maps[['Afrocarpus usambarensis']]$scenario
    domains <- c('narrow', 'broad')
    paired <- all(vapply(domains, function(domain) {
      moderate <- afrocarpus[[paste0(domain, '_moderate')]]
      stronger <- afrocarpus[[paste0(domain, '_stronger')]]
      !is.null(moderate) && !is.null(stronger) &&
        identical(terra::values(moderate), terra::values(stronger))
    }, logical(1)))
    treatments <- c(baseline = 'Unthinned', moderate = '8 km thinning', stronger = '16 km thinning')
    treatments <- treatments[names(treatments) %in% sub('^[^_]+_', '', names(afrocarpus))]
    if (paired) treatments <- c(baseline = 'Unthinned', moderate = '8 and 16 km thinning')

    columns <- length(treatments)
    row_heights <- c(.14, .14, 1, .14, .14, 1)
    plot_height <- par('din')[2] - sum(par('omi')[c(1, 3)])
    extra_margin <- (3.1 - .55) * par('csi') * par('mex')
    extra_bottom_height <- extra_margin * (sum(row_heights) + 1) / (plot_height - extra_margin)

    layout(rbind(
      c(0, rep(1, columns)),
      c(0, 2, 3, rep(0, columns - 2)),
      c(0, 4, 5, rep(0, columns - 2)),
      c(0, rep(6, columns)),
      c(0, seq_len(columns) + 6),
      matrix(seq.int(7 + columns, length.out = 2 * (columns + 1)), nrow = 2, byrow = TRUE)
    ), widths = c(.13, rep(1, columns)), heights = c(row_heights, 1 + extra_bottom_height))

    heading <- function(label, strip = NULL, vertical = FALSE, left_margin = 3.15,
                        thin = FALSE, bottom_margin = .55, label_cex = 18/12) {
      par(mar = if (vertical) c(bottom_margin, 0, .25, 0) else c(0, left_margin, 0, .45))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
      if (!is.null(strip)) {
        x <- if (vertical) c(.25, 1) else c(0, 1)
        y <- if (vertical) c(0, 1) else if (thin) c(.16, .84) else c(.08, .92)
        rect(x[1], y[1], x[2], y[2], col = strip, border = 'grey50', xpd = TRUE)
      }
      text(if (vertical) .625 else .5, .5, label, cex = label_cex, srt = if (vertical) 90 else 0)
    }
    panel <- 0L
    frames <- list()
    extent <- as.vector(terra::ext(afrocarpus[['broad_baseline']]))
    rank_map <- function(raster, show_latitude) {
      bottom_row <- panel >= 2L + columns
      par(mar = c(if (bottom_row) 3.1 else .55, if (show_latitude) 3.15 else .8, .25, .45))
      plot.new()
      reference_width <- (par('din')[1] - sum(par('omi')[c(2, 4)])) /
        (.13 + columns) - (3.15 + .45) * par('csi') * par('mex')
      x_half <- (diff(extent[1:2]) + .30) / 2 * par('pin')[1] / reference_width
      plot.window(xlim = mean(extent[1:2]) + c(-x_half, x_half),
        ylim = extent[3:4] + c(-.15, .15), asp = 1, xaxs = 'i', yaxs = 'i')
      terra::plot(raster, col = suit, range = c(0, 1), legend = FALSE, add = TRUE)
      if (bottom_row) {
        axis(1, at = seq(34, 42, 2), cex.axis = 17/12, mgp = c(2, .65, 0), col = frame_colour)
      }
      if (show_latitude) {
        axis(2, at = seq(-4, 4, 2), las = 1, cex.axis = 17/12, mgp = c(2, .65, 0), col = frame_colour)
      }
      box(col = frame_colour)
      panel <<- panel + 1L
      bounds <- par('usr')
      frames[[panel]] <<- c(grconvertX(bounds[1:2], 'user', 'ndc'),
        grconvertY(bounds[3:4], 'user', 'ndc'))
      text(bounds[1] + .012 * diff(bounds[1:2]), bounds[4] - .012 * diff(bounds[3:4]),
        paste0('(', letters[panel], ')'), adj = c(0, 1), font = 2, cex = 18/12)
    }
    heading(expression(bolditalic(Coffea~fadenii)), strip = '#E5E5E5', thin = TRUE, label_cex = 16/12)
    for (domain in domains) heading(paste(tools::toTitleCase(domain), 'calibration'),
      strip = '#E5E5E5', left_margin = if (domain == 'narrow') 3.15 else .8,
      thin = TRUE)
    for (domain in domains) rank_map(coffea[[paste0(domain, '_baseline')]], domain == 'narrow')
    heading(expression(bolditalic(Afrocarpus~usambarensis)), strip = '#E5E5E5', thin = TRUE, label_cex = 16/12)
    for (i in seq_along(treatments)) heading(treatments[i],
      strip = '#E5E5E5', left_margin = if (i == 1L) 3.15 else .8, thin = TRUE)
    for (domain in domains) {
      heading(paste(tools::toTitleCase(domain), 'calibration'), strip = '#E5E5E5', vertical = TRUE, bottom_margin = if (domain == 'broad') 3.1 else .55)
      for (treatment in names(treatments)) {
        rank_map(afrocarpus[[paste(domain, treatment, sep = '_')]], treatment == names(treatments)[1])
      }
    }

    first_bottom <- 3L + columns
    left_edge <- frames[[first_bottom]][1]
    right_edge <- frames[[length(frames)]][2]
    device_width <- par('din')[1]
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
    breaks <- seq(0, 1, length.out = length(suit) + 1)
    legend_x <- grconvertX(right_edge + c(.12, .25) / device_width, 'ndc', 'user')
    legend_rows <- list(frames[[2]][3:4],
      c(frames[[length(frames)]][3], frames[[2 + columns]][4]))
    for (edges in legend_rows) {
      y <- grconvertY(edges[1] + breaks * diff(edges), 'ndc', 'user')
      rect(legend_x[1], head(y, -1), legend_x[2], tail(y, -1),
        col = suit, border = NA, xpd = NA)
      rect(legend_x[1], y[1], legend_x[2], tail(y, 1), border = 'black', xpd = NA)
      tick_y <- seq(y[1], tail(y, 1), length.out = 5)
      text(grconvertX(right_edge + .30 / device_width, 'ndc', 'user'),
        tick_y, paste0(seq(0, 100, 25), '%'), adj = c(0, .5), cex = 17/12, xpd = NA)
    }
    longitude_y <- frames[[first_bottom]][3] - .42 / par('din')[2]
    text(grconvertX(mean(c(left_edge, right_edge)), 'ndc', 'user'),
      grconvertY(longitude_y, 'ndc', 'user'), 'Longitude (\u00b0)',
      font = 2, cex = 18/12, xpd = NA)
    for (y in c(mean(frames[[1]][3:4]),
      mean(c(frames[[3]][4], frames[[first_bottom]][3])))) {
      latitude_strip_edge <- par('omi')[2] + .13 / (.13 + columns) * (device_width - sum(par('omi')[c(2, 4)]))
      latitude_half_width <- strheight('Latitude (\u00b0)', units = 'inches',
        font = 2, cex = 18/12) / 2
      latitude_x <- (latitude_strip_edge + 2/25.4 + latitude_half_width) / par('din')[1]
      text(grconvertX(latitude_x, 'ndc', 'user'), grconvertY(y, 'ndc', 'user'),
        'Latitude (\u00b0)', srt = 90, font = 2, cex = 18/12, xpd = NA)
    }
  }
  draw_differences <- function() {
    par(mar = rep(0, 4), oma = rep(0, 4), bg = 'white', xpd = NA)
    colours <- c('#E6E6E6', '#D55E00', '#0072B2', '#00856A')
    device <- par('din')
    half_gap <- 1 / (25.4 * device[1])
    width <- 6.19344 / device[1]
    centre <- 1 - width - half_gap - .10 / device[1]
    vertical <- c(.24, .92)
    extent <- Reduce(terra::union, lapply(maps, function(m) terra::ext(m$difference)))
    e <- as.vector(extent)
    for (i in seq_along(maps)) {
      m <- maps[[i]]
      horizontal <- if (i == 1L) centre - half_gap - c(width, 0) else centre + half_gap + c(0, width)
      par(fig = c(0, 1, 0, 1), plt = c(horizontal, vertical), new = i != 1L)
      plot.new()
      plot.window(xlim = range(e[1:2] + c(-.15, .15), 42.15),
        ylim = e[3:4] + c(-.15, .15), asp = 1, xaxs = 'i', yaxs = 'i')
      draw_raster(m$difference, colours, seq(-.5, 3.5))
      hatch(m$novel, density = 8, col = '#777777')
      map_axes(i, TRUE)
      map_heading(letters[i], titles[[i]][[1]],
        paste(scenario_label(m$difference_names), collapse = ' vs. '), species_below = TRUE)
    }
    par(fig = c(0, 1, 0, 1), mar = rep(0, 4), new = TRUE)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
    text(centre, .155, 'Longitude (\u00b0)', font = 2, cex = 17/12)
    text(.12 / device[1], mean(vertical), 'Latitude (\u00b0)', srt = 90, font = 2, cex = 17/12)
    labels <- c('First scenario', 'Second scenario', 'Both', 'Neither')
    key_width <- .23 / device[1]
    key_height <- .20 / device[2]
    label_gap <- .09 / device[1]
    entry_widths <- key_width + label_gap + strwidth(labels, cex = 15/12)
    entry_gap <- .25 / device[1]
    left <- (centre - (sum(entry_widths) + 3 * entry_gap)/2) +
      c(0, cumsum(head(entry_widths, -1) + entry_gap))
    rect(left, .095 - key_height/2, left + key_width, .095 + key_height/2,
      col = colours[c(3, 2, 4, 1)], border = 'black')
    text(left + key_width + label_gap, .095, labels, adj = c(0, .5), cex = 15/12)
    rect(left[1], .032 - key_height/2, left[1] + key_width, .032 + key_height/2,
      density = 8, angle = 45, col = '#777777', border = 'black')
    text(left[1] + key_width + label_gap, .032, 'Outside training ranges',
      adj = c(0, .5), cex = 15/12)
  }
  figures <- list(Figure_7_Maxent_predictions_of_climatic_suitability_in_Kenya = draw_main,
    Figure_S9_Occurrence_cells_and_calibration_areas = draw_occurrences,
    Figure_S10_Percentile_ranks_of_predicted_climatic_suitability = draw_scenarios,
    Figure_S11_Agreement_between_selected_Maxent_scenarios = draw_differences)
  for (n in names(figures))for (ext in c('pdf', 'png')) {
    path <- rq_file(out, paste0(n, '.', ext))
    width <- if (n == 'Figure_7_Maxent_predictions_of_climatic_suitability_in_Kenya')15.18 else if (n == 'Figure_S10_Percentile_ranks_of_predicted_climatic_suitability')10 else if (n == 'Figure_S11_Agreement_between_selected_Maxent_scenarios')13.30 else 11
    height <- if (n == 'Figure_7_Maxent_predictions_of_climatic_suitability_in_Kenya')14 else if (n == 'Figure_S11_Agreement_between_selected_Maxent_scenarios')9 else if (n == 'Figure_S10_Percentile_ranks_of_predicted_climatic_suitability')10.2 else 9
    if (n %in% c('Figure_7_Maxent_predictions_of_climatic_suitability_in_Kenya', 'Figure_S9_Occurrence_cells_and_calibration_areas', 'Figure_S10_Percentile_ranks_of_predicted_climatic_suitability', 'Figure_S11_Agreement_between_selected_Maxent_scenarios')) {
      if (ext == 'pdf') cairo_pdf(path, width = width, height = height, family = 'Arial')
      else png(path, width = width * 200, height = height * 200, res = 200,
        type = 'cairo', family = 'Arial')
    } else if (ext == 'pdf')pdf(path, width = width, height = height, compress = dir.exists(tempdir()))
    else png(path, width = width * 200, height = height * 200, res = 200)
    tryCatch(figures[[n]](), finally = dev.off())
  }
}

# ============================================================
# 18. Hs occurrence maps
# ============================================================

# Plot source evidence while the inputs for crossed Maxent remain under review.
hs_occurrence_maps <- function(root) {
  .libPaths(c(file.path(root, 'analysis/scripts/.r-library'), .libPaths()))
  if (!requireNamespace('terra', quietly = TRUE))stop('Required package unavailable: terra')
  out <- rq_dirs(root, 5)
  read <- function(x)hs_read(file.path(root, x))
  ledger <- read('analysis/Datasets/derived/worksheet_focal_count_ledger.csv')
  coords <- read('analysis/Datasets/derived/worksheet_coordinate_transcription.csv')
  refs <- read('analysis/Datasets/derived/study_site_reference_coordinates.csv')
  stopifnot(!anyDuplicated(coords$plot_id))
  species <- c('Coffea fadenii', 'Afrocarpus usambarensis')
  county <- terra::vect(file.path(root, 'analysis/Datasets/mapping/boundaries/ken_admbnda_adm1_iebc_20191031.shp'))
  county <- terra::project(county, 'EPSG:4326')
  county_names <- terra::values(county)
  is_taita <- apply(county_names, 1, function(x)any(grepl('taita.*taveta', x, ignore.case = TRUE)))
  stopifnot(sum(is_taita) == 1)
  taita <- county[is_taita,]
  field <- gbif <- list()
  for (sp in species) {
    z <- ledger[which(ledger$species == tolower(sp) & hs_bool(ledger$included) &
        hs_bool(ledger$stage_table_available) & is.finite(ledger$count) & ledger$count>0),]
    f <- unique(z[c('record_id', 'forest', 'workbook', 'sheet')])
    m <- match(f$record_id, coords$plot_id)
    f$species <- sp
    f$lon <- coords$source_lon[m]
    f$lat <- coords$source_lat[m]
    f$coordinate_status <- ifelse(hs_coord_ok(f$lon, f$lat),
      'Mapped worksheet coordinate; accuracy and original datum unverified', 'No usable coordinate')
    field[[sp]] <- f
    g <- read(file.path('analysis/Datasets/mapping', paste0(tolower(sp), '_GBIF.csv')))
    names(g)[names(g) == 'decimalLongitude'] <- 'lon'
    names(g)[names(g) == 'decimalLatitude'] <- 'lat'
    names(g)[names(g) == 'coordinateUncertaintyInMeters'] <- 'uncertainty_m'
    g$map_status <- ifelse(hs_coord_ok(g$lon, g$lat) & g$species %in% sp &
      g$occurrenceStatus %in% 'PRESENT' & g$countryCode %in% c('KE', 'TZ'),
      'Mapped source record; taxonomy, wild status and coordinate accuracy pending review',
      'Outside displayed source scope or missing usable presence coordinates')
    gbif[[sp]] <- g
  }
  field_all <- do.call(rbind, field)
  gbif_all <- do.call(rbind, gbif)
  write.csv(field_all, rq_file(out, 'mapped_field_occurrence_evidence.csv'), row.names = FALSE, na = 'NA')
  write.csv(gbif_all, rq_file(out, 'GBIF_occurrence_map_review.csv'), row.names = FALSE, na = 'NA')
  local_x <- range(c(field_all$lon, refs$long), na.rm = TRUE) + c( - .024, .024)
  local_y <- range(c(field_all$lat, refs$lat), na.rm = TRUE) + c( - .023, .023)
  panels <- list(expression(italic(Coffea ~ fadenii)), expression(italic(Afrocarpus ~ usambarensis)))
  draw <- function() {
    par(mfrow = c(2, 2), mar = c(3.7, 4.3, 3, 1), oma = c(3, 0, 1, 0))
    for (i in seq_along(species)) {
      sp <- species[i]
      f <- field[[sp]]
      f <- f[hs_coord_ok(f$lon, f$lat),]
      g <- gbif[[sp]]
      g <- g[startsWith(g$map_status, 'Mapped'),]
      for (local in c(FALSE, TRUE)) {
        x <- if (local)local_x else range(c(g$lon, local_x)) + c( - .35, .35)
        y <- if (local)local_y else range(c(g$lat, local_y)) + c( - .35, .35)
        plot(NA, xlim = x, ylim = y, asp = 1, las = 1, xlab = '',
          ylab = 'Latitude (degrees; negative = south)', main = panels[[i]])
        terra::plot(county, add = TRUE, col = '#F4F4F4', border = '#C4C4C4', lwd = .5)
        terra::plot(taita, add = TRUE, col = '#E7EFE5', border = '#527451', lwd = 1)
        points(g$lon, g$lat, pch = 1, col = '#8C4A83', cex = .85)
        points(f$lon, f$lat, pch = 17, col = '#0072B2', cex = .9)
        if (local) {
          points(refs$long, refs$lat, pch = 4, lwd = 1.5)
          text(refs$long, refs$lat, sub(' forest$', '', refs$study_site), pos = 1, cex = .8)
        } else rect(local_x[1], local_y[1], local_x[2], local_y[2], border = '#4C754A', lwd = 1.4, lty = 2)
        box()
        mtext('Longitude (degrees E)', side = 1, line = 2.15)
        mtext(if (local)'Taita Hills detail'else 'Kenya and Tanzania records', side = 3, line = .25, cex = .8)
      }
    }
    mtext('Blue triangles: worksheet detections   |   Purple circles: archived GBIF   |   Crosses: approximate forest references',
      side = 1, outer = TRUE, line = .8, cex = .8)
  }
  for (ext in c('pdf', 'png')) {
    path <- rq_file(out, paste0('Taita_focal_species_occurrence_evidence.', ext))
    if (ext == 'pdf')pdf(path, width = 11, height = 9, compress = dir.exists(tempdir()))
    else png(path, width = 2200, height = 1800, res = 200)
    tryCatch(draw(), finally = dev.off())
  }
  writeLines(c('Occurrence evidence map for both focal species. The map displays primary worksheet records with positive listed evidence and archived Kenya/Tanzania GBIF records with usable presence coordinates. Duplicate or nearby source points can overlap.',
      'The green fill shows Taita-Taveta County and the dashed rectangle locates the local panel. Crosses are approximate forest reference points. The local panel is a display window.',
      'GBIF metadata and worksheet coordinates remain subject to review. These records have not been certified as Maxent inputs. This figure depicts source evidence only.'), rq_file(out, 'OCCURRENCE_MAP_CAPTION.txt'))
  invisible(out)
}

# ============================================================
# 19. Hs prepare taita exploratory
# ============================================================

hs_prepare_taita_exploratory <- function(root, archive) {
  .libPaths(c(file.path(root, 'analysis/scripts/.r-library'), .libPaths()))
  if (!requireNamespace('terra', quietly = TRUE))stop('Required package unavailable: terra')
  if (!file.exists(archive))stop('Provide the WorldClim 2.1 2.5-minute BIO archive through TAITA_WORLDCLIM_ZIP.')
  input <- file.path(root, 'analysis/Datasets/habitat_suitability/exploratory')
  dir.create(input, recursive = TRUE, showWarnings = FALSE)
  read <- function(path)hs_read(file.path(root, path))
  save <- function(x, name)write.csv(x, file.path(input, name), row.names = FALSE, na = 'NA')
  sp <- c('Coffea fadenii', 'Afrocarpus usambarensis')
  ledger <- read('analysis/Datasets/derived/worksheet_focal_count_ledger.csv')
  coordinates <- read('analysis/Datasets/derived/worksheet_coordinate_transcription.csv')
  refs <- read('analysis/Datasets/derived/study_site_reference_coordinates.csv')
  save(data.frame(label = sub(' forest$', '', refs$study_site), lon = refs$long, lat = refs$lat),
    'forest_reference_points.csv')
  records <- lapply(sp, function(species) {
      z <- ledger[which(ledger$species == tolower(species) & hs_bool(ledger$included) &
          hs_bool(ledger$stage_table_available) & is.finite(ledger$count) & ledger$count>0),]
      id <- unique(z$record_id)
      at <- match(id, coordinates$plot_id)
      field <- data.frame(species, source = 'Field', source_id = id, lon = coordinates$source_lon[at],
        lat = coordinates$source_lat[at], uncertainty_m = NA_real_, collection_year = NA_integer_,
        duplicate_of = NA_character_, screen_pass = hs_coord_ok(coordinates$source_lon[at], coordinates$source_lat[at]),
        source_detail = 'Listed positive counts in primary forest-associated worksheet')
      g <- read(file.path('analysis/Datasets/mapping', paste0(tolower(species), '_GBIF.csv')))
      bad <- grepl('COORDINATE_INVALID|COORDINATE_OUT_OF_RANGE|COUNTRY_COORDINATE_MISMATCH|ZERO_COORDINATE',
        g$issue, ignore.case = TRUE)
      managed <- grepl('cultivat|garden|introduced|plantation', paste(g$establishmentMeans, g$locality), ignore.case = TRUE)
      dated <- is.na(g$year) | g$year >= 1900
      screen <- hs_coord_ok(g$decimalLongitude, g$decimalLatitude) & g$species %in% species &
      g$countryCode %in% c('KE', 'TZ') & g$occurrenceStatus %in% 'PRESENT' & !bad & !managed & dated
      occurrence_key <- trimws(as.character(g$occurrenceID))
      key_rows <- which(screen & !is.na(occurrence_key) & nzchar(occurrence_key))
      first <- match(occurrence_key[key_rows], occurrence_key[key_rows])
      duplicate_of <- rep(NA_character_, nrow(g))
      repeated <- which(first != seq_along(first))
      duplicate_of[key_rows[repeated]] <- as.character(g$gbifID[key_rows[first[repeated]]])
      gbif <- data.frame(species, source = 'GBIF', source_id = as.character(g$gbifID),
        lon = g$decimalLongitude, lat = g$decimalLatitude, uncertainty_m = g$coordinateUncertaintyInMeters,
        collection_year = g$year, duplicate_of, screen_pass = screen,
        source_detail = paste(g$countryCode, g$basisOfRecord, g$issue, sep = ' | '))
      rbind(field, gbif)
    })
  r <- do.call(rbind, records)
  r$audit_id <- paste(r$species, r$source, r$source_id, sep = '|')
  r$taxonomy_verified <- r$coordinate_verified <- r$wild_verified <- FALSE
  save(r, 'screened_occurrences.csv')
  members <- paste0('wc2.1_2.5m_bio_', c(1, 12, 15), '.tif')
  td <- tempfile('taita-worldclim-');
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  unzip(archive, files = members, exdir = td)
  kept <- r[r$screen_pass & is.na(r$duplicate_of),]
  stopifnot(nrow(kept)>0)
  boundary_path <- file.path(root, 'analysis/Datasets/mapping/boundaries/ken_admbnda_adm0_iebc_20191031.shp')
  kenya <- terra::project(terra::vect(boundary_path), 'EPSG:4326')
  bounds <- as.vector(terra::ext(kenya))
  x <- range(range(kept$lon) + c( - 2.2, 2.2), bounds[1:2])
  y <- range(range(kept$lat) + c( - 2.2, 2.2), bounds[3:4])
  env <- terra::crop(terra::rast(file.path(td, members)), terra::ext(x[1], x[2], y[1], y[2]))
  names(env) <- c('bio1', 'bio12', 'bio15')
  terra::writeRaster(env, file.path(input, 'baseline.tif'), overwrite = TRUE)
  xy <- terra::xyFromCell(env, seq_len(terra::ncell(env)))
  usable <- complete.cases(terra::values(env))
  mask <- function(included, name) {
    z <- env[[1]]
    terra::values(z) <- ifelse(included & usable, 1, NA_real_)
    terra::writeRaster(z, file.path(input, name), overwrite = TRUE)
  }
  for (species in sp) {
    p <- kept[kept$species == species,]
    stopifnot(nrow(p)>0)
    distance <- rep(Inf, nrow(xy))
    for (i in seq_len(nrow(p)))distance <- pmin(distance, mx_distance(xy[, 1], xy[, 2], p$lon[i], p$lat[i]))
    stem <- gsub(' ', '_', species)
    mask(distance <= 75000, paste0(stem, '_narrow.tif'))
    mask(distance <= 175000, paste0(stem, '_broad.tif'))
  }
  national <- terra::rasterize(kenya, env[[1]], field = 1, background = NA, touches = FALSE)
  mask(!is.na(terra::values(national, mat = FALSE)), 'prediction_mask.tif')
  save(data.frame(species = sp, thinning_m = 8000, stronger_thinning_m = 16000,
      thinning_rationale = '8 km is above the ~4.6 km predictor cell width and addresses adjacent-cell clustering.',
      stronger_thinning_rationale = '16 km tests greater separation; omit if any repeat has fewer than six cells.',
      narrow_rationale = '75 km buffers around screened localities, limited to valid climate cells. Exploratory proximity assumption.',
      broad_rationale = '175 km buffers around the same localities. Exploratory wider background assumption.',
      prediction_rationale = 'Kenya national boundary (IEBC 2019). Include predictor cells whose centres lie within Kenya and have all three climate values.'),
    'crossed_design.csv')
  writeLines(c('WorldClim 2.1 1970-2000 normals, 2.5 arc-minutes. Source: https://worldclim.org/data/worldclim21.html',
      paste('Archive MD5:', tools::md5sum(archive)),
      'Prediction mask uses the Kenya IEBC national boundary (2019), retaining climate-valid cell centres inside the country.',
      'GBIF screen keeps usable Kenya/Tanzania presence coordinates and expected accepted name, excludes specified coordinate issues and indications of cultivation or introduction, and retains unknown dates and uncertainty. Known dates before 1900 are excluded. Occurrence IDs and environmental cells are deduplicated.',
      'Field points use listed-positive primary worksheets. Original datum, wild context and positional uncertainty are not independently verified. Field and GBIF verification flags remain false.',
      'Calibration buffers are scenario proxies, not validated accessible areas. Predictions are regional climate rankings at native resolution.'), file.path(input, 'PROVENANCE.txt'))
  input
}

# ============================================================
# 20. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep('^--file=', commandArgs(), value = TRUE)
  start <- if (length(arg))dirname(normalizePath(sub('^--file=', '', arg[1]), mustWork = TRUE))else getwd()
  root <- normalizePath(start, winslash = '/', mustWork = TRUE)
  while (!file.exists(file.path(root, 'analysis/scripts/05_climatic_suitability.R'))) {
    parent <- dirname(root)
    if (identical(parent, root))stop('Cannot locate the Taita project from ', start, call. = FALSE)
    root <- parent
  }
  input_root <- file.path(root, 'analysis/Datasets/habitat_suitability')
  screened_dir <- file.path(input_root, 'exploratory')
  mode <- Sys.getenv('TAITA_HABITAT_MODE')
  if (!nzchar(mode)) {
    mode <- 'maxent_exploratory'
    if (!nzchar(Sys.getenv('TAITA_HABITAT_INPUT_DIR')) &&
      !all(file.exists(file.path(input_root, hs_crossed_inputs())))) {
      if (all(file.exists(file.path(screened_dir, hs_crossed_inputs(screened = TRUE)))) ||
        nzchar(Sys.getenv('TAITA_WORLDCLIM_ZIP')))mode <- 'taita_screened'
    }
  }
  message('Question 5 mode: ', mode)
  if (mode == 'maxent_exploratory') {
    run_maxent_exploratory(root)
  } else if (mode == 'taita_screened') {
    archive <- Sys.getenv('TAITA_WORLDCLIM_ZIP')
    design_file <- file.path(screened_dir, 'crossed_design.csv')
    old_extent <- file.exists(design_file) && any(grepl('Fixed Taita Hills window', readLines(design_file, warn = FALSE)))
    if (!nzchar(archive) && old_extent)
    archive <- file.path(input_root, 'source/taita_wc21_selected_2.5m.zip')
    if (nzchar(archive))screened_dir <- hs_prepare_taita_exploratory(root, archive)
    run_maxent_exploratory(root, screened_dir, screened = TRUE)
  } else if (mode == 'occurrence_maps') {
    hs_occurrence_maps(root)
  } else if (mode == 'cached_maps') {
    run_habitat_suitability(root)
  } else stop('Unknown habitat mode: ', mode, call. = FALSE)
  source(file.path(file.path(root, 'analysis'), 'scripts', 'combine_output_tables.R'), local = TRUE)
  combine_output_tables(file.path(file.path(root, 'analysis'), 'outputs', 'research question 5'), '05_climatic_suitability.R')
}
