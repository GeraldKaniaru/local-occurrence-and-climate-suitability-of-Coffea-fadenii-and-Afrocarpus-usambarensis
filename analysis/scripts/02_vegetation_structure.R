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
# 1. Vs project root
# ============================================================

# Research question 2: How does vegetation structure differ between the surveyed neighbourhoods of the two fragments, and which woody taxa have the highest IVI in each fragment?
# Descriptive vegetation structure in surveyed Taita neighbourhoods.

vs_project_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = '/', mustWork = TRUE)
  repeat {if (file.exists(file.path(p, 'analysis/scripts/02_vegetation_structure.R')))return(p)
    parent <- dirname(p);
    if (identical(parent, p))break;
    p <- parent}
  stop('Set the working directory to the Taita project or a folder within it.', call. = FALSE)
}

# ============================================================
# 2. Vs temp ready
# ============================================================

vs_temp_ready <- function() {
  tmp <- tempdir(check = TRUE);
  if (!dir.exists(tmp))dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  probe <- tempfile('vegetation_write_check_', tmpdir = tmp);
  ok <- suppressWarnings(file.create(probe))
  if (isTRUE(ok))unlink(probe)else stop('R temporary directory is unavailable. Restart R and Source again.', call. = FALSE)
  invisible(tmp)
}

# ============================================================
# 3. Vs table 2a
# ============================================================

vs_table_2a <- function(summary) {
  metrics <- c('listed_stems', 'median_dbh_cm', 'quadratic_mean_dbh_cm', 'median_height_m',
    'observed_basal_area_m2', 'canopy_percent', 'shrub_percent', 'herb_percent', 'litter_percent', 'bare_percent')
  measures <- c('Listed stems per record', 'Median DBH per record (cm)',
    'Quadratic mean DBH per record (cm)', 'Median height per record (m)',
    'Basal area per record (m\u00b2)', 'Canopy cover (%)', 'Shrub cover (%)',
    'Herbaceous cover (%)', 'Litter cover (%)', 'Bare soil and rock cover (%)')
  digits <- c(1L, 2L, 2L, 2L, 2L, 1L, 0L, 0L, 0L, 0L)
  table <- data.frame(Measure = measures)
  for (forest in c('Ngangao', 'Mbololo')) {
    d <- summary[which(summary$forest == forest),]
    d <- d[match(metrics, d$metric),]
    format_value <- function(x) sprintf('%.*f', digits, floor(x * 10^digits + .5) / 10^digits)
    table[[forest]] <- ifelse(complete.cases(d[c('median', 'q1', 'q3')]),
      sprintf('%s [%s, %s]', format_value(d$median), format_value(d$q1), format_value(d$q3)), NA_character_)
  }
  table
}

# ============================================================
# 4. Export vegetation reporting
# ============================================================

export_vegetation_reporting <- function(out) {
  read <- function(n)read.csv(rq_file(out, paste0(n, '.csv')), check.names = FALSE)
  forest <- read('forest_record_summary');
  forest <- forest[forest$cohort == 'primary',]
  cv <- read('forest_cover_summary')
  metrics <- c(listed_woody_stems = 'listed_stems', dbh_median_cm = 'median_dbh_cm',
    height_median_m = 'median_height_m', summed_basal_area_m2 = 'observed_basal_area_m2')
  matched <- forest$metric %in% names(metrics)
  forest$metric[matched] <- unname(metrics[forest$metric[matched]])
  cover_metrics <- c('Upper canopy cover' = 'canopy_percent', 'Shrub layer cover' = 'shrub_percent',
    'Herbaceous layer cover' = 'herb_percent', 'Litter cover' = 'litter_percent', 'Bare ground and rocks' = 'bare_percent')
  summary <- rbind(forest[c('forest', 'metric', 'median', 'q1', 'q3')],
    data.frame(forest = cv$forest, metric = unname(cover_metrics[cv$attribute]),
      median = cv$median_percent, q1 = cv$q1_percent, q3 = cv$q3_percent))
  tables <- list(Table_S4_Vegetation_structure_and_cover = vs_table_2a(summary),
    Table_S4_record_level_vegetation_structure = read('vegetation_structure_by_record'),
    Table_S5_measurement_completeness = read('measurement_completeness'),
    Table_S5_Importance_Value_Index_among_listed_larger_stems = vs_table_s3(subset(read('importance_value_index'), cohort == 'primary')))
  for (n in names(tables))write.csv(tables[[n]], rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  figures <- c('Figure_2_recorded_vegetation_structure', 'Figure_S2_vegetation_layers_and_ground_cover',
    'Figure_S3_pooled_diameter_and_height_distributions', 'Figure_S4_importance_value_index_components')
  write.csv(data.frame(label = c('Figure 2', 'Figure S2', 'Figure S3', 'Figure S4', 'Table 2', 'Table S4', 'Table S5', 'Table S3'),
      file = rq_relative(c(paste0(figures, '.pdf'), paste0(names(tables), '.csv')))), rq_file(out, 'manuscript_reporting_register.csv'), row.names = FALSE)
  writeLines(c('# Vegetation-structure reporting', '',
      'Figure 2: Record-level stem abundance, diameter, quadratic mean diameter, summed individual basal area and height. Points represent records; forest inference is limited to surveyed neighbourhoods.',
      'Figure S2: Separate vegetation-layer and ground-surface cover attributes.',
      'Figure S3: Pooled diameter and height distributions with measurement denominators; records with more listed stems contribute more heavily.',
      'Figure S4: Top-ranked taxa by descriptive IVI, showing relative stem abundance, relative record frequency and relative basal area.',
      'Table 2: Record-level structure and cover medians and quartiles by forest.',
      'Table 3: Leading taxa and the three IVI components by forest.',
      'Table S4: Record-level structural metrics and eligibility.',
      'Table S5: Measurement completeness, missingness and invalid values.',
      'Table S3: Full IVI rankings and component denominators.',
      'Source files remain in existing relevant analysis folders; no extra input or helper folders are required.'), rq_file(out, 'REPORTING_GUIDE.md'))
}

# ============================================================
# 5. Vs size summary
# ============================================================

vs_size_summary <- function(x) {
  valid <- x[is.finite(x) & x>0]
  if (!length(valid)) return(c(n = 0, median = NA_real_, q1 = NA_real_, q3 = NA_real_,
      iqr = NA_real_, minimum = NA_real_, maximum = NA_real_))
  qs <- unname(quantile(valid, c(.25, .5, .75), type = 7))
  c(n = length(valid), median = qs[2], q1 = qs[1], q3 = qs[3], iqr = qs[3] - qs[1],
    minimum = min(valid), maximum = max(valid))
}

# ============================================================
# 6. Vs basal area
# ============================================================

vs_basal_area <- function(dbh_cm) {
  ifelse(is.finite(dbh_cm) & dbh_cm>0, pi * (dbh_cm / 200)^2, NA_real_)
}

# ============================================================
# 7. Vs table s3
# ============================================================

vs_table_s3 <- function(ivi) {
  ivi <- ivi[order(ivi$forest, ivi$rank),]
  decimal <- function(x) ifelse(is.na(x), NA_character_, sprintf('%.2f', x))
  data.frame(
    'Forest' = ivi$forest,
    'Rank' = ivi$rank,
    'Taxon' = sub('^Ocotea usambarensis$', 'Kuloa usambarensis',
      paste0(toupper(substr(ivi$taxon, 1, 1)), substring(ivi$taxon, 2)), ignore.case = TRUE),
    'Stems' = ivi$listed_stems,
    'Relative abundance (%)' = decimal(ivi$relative_stem_abundance_percent),
    'Relative frequency (%)' = decimal(ivi$relative_record_frequency_percent),
    'Relative basal area (%)' = decimal(ivi$relative_dominance_percent),
    'IVI points' = decimal(ivi$IVI_points),
    check.names = FALSE)
}

# ============================================================
# 8. Vs ivi
# ============================================================

vs_ivi <- function(stems, record_ids, forest, cohort) {
  z <- stems[stems$plot_id %in% record_ids & stems$woody_eligible,, drop = FALSE]
  if (!nrow(z))return(data.frame())
  taxa <- sort(unique(z$analysis_name))
  out <- do.call(rbind, lapply(taxa, function(taxon) {
        s <- z[z$analysis_name == taxon,, drop = FALSE]
        data.frame(cohort, forest, taxon, listed_stems = nrow(s),
          records_with_taxon = length(unique(s$plot_id)), survey_records = length(record_ids),
          dbh_available = sum(is.finite(s$basal_area_m2)),
          dbh_missing_or_invalid = sum(!is.finite(s$basal_area_m2)),
          observed_basal_area_m2 = sum(s$basal_area_m2, na.rm = TRUE))
      }))
  out$forest_listed_stem_total <- sum(out$listed_stems)
  out$summed_taxon_record_occurrences <- sum(out$records_with_taxon)
  out$forest_observed_basal_area_m2 <- sum(out$observed_basal_area_m2)
  out$basal_area_complete <- all(out$dbh_missing_or_invalid == 0)
  out$record_frequency_percent <- 100 * out$records_with_taxon / out$survey_records
  out$relative_stem_abundance_percent <- 100 * out$listed_stems / out$forest_listed_stem_total
  out$relative_record_frequency_percent <- 100 * out$records_with_taxon / out$summed_taxon_record_occurrences
  out$relative_observed_basal_area_percent <- if (out$forest_observed_basal_area_m2[1]>0)
  100 * out$observed_basal_area_m2 / out$forest_observed_basal_area_m2 else NA_real_
  # Missing size data must not silently turn a partial BA inventory into a full
  # dominance component. Keep observed shares for audit; withhold complete IVI.
  complete <- out$basal_area_complete[1] && out$forest_observed_basal_area_m2[1]>0
  out$relative_dominance_percent <- if (complete)out$relative_observed_basal_area_percent else NA_real_
  out$IVI_points <- out$relative_stem_abundance_percent +
  out$relative_record_frequency_percent + out$relative_dominance_percent
  out$IVI_percent <- out$IVI_points / 3
  out$rank <- if (complete)rank( - out$IVI_points, ties.method = 'min')else NA_integer_
  out$status <- if (complete)'Descriptive IVI from listed woody stems; field units and sampling unverified'else
  'IVI unavailable: incomplete or nonpositive forest basal-area denominator'
  out[order(out$rank, out$taxon, na.last = TRUE),, drop = FALSE]
}

# ============================================================
# 9. Vs ranked ivi plot
# ============================================================

vs_ranked_ivi_plot <- function(ivi) {
  forests <- c('Ngangao', 'Mbololo')
  components <- c('Relative abundance', 'Relative frequency', 'Relative basal area')
  columns <- c('relative_stem_abundance_percent', 'relative_record_frequency_percent',
    'relative_dominance_percent')
  top <- lapply(forests, function(f) {
    z <- ivi[which(ivi$forest == f & is.finite(ivi$rank)), ]
    head(z[order(z$rank, z$taxon), ], 10)
  })
  limit <- 50
  palette <- c('#0072B2', '#E69F00', '#009E73')
  panels <- lapply(seq_along(top), function(i) {
    z <- top[[i]]
    names <- paste0(toupper(substr(z$taxon, 1, 1)), substring(z$taxon, 2))
    names <- sub('^Ocotea usambarensis$', 'Kuloa usambarensis', names, ignore.case = TRUE)
    values <- data.frame(taxon = factor(rep(names, 3), levels = rev(unique(names))),
      component = factor(rep(components, each = nrow(z)), levels = components),
      value = unlist(z[columns], use.names = FALSE))
    values$endpoint <- rep(z$IVI_points, 3)
    values$label_y <- as.integer(values$taxon) +
      c(.16, -.16, .16)[as.integer(values$component)]
    values$separator_x <- values$endpoint + 7.5
    values$label_x <- values$separator_x +
      c(-.75, -.75, .75)[as.integer(values$component)]
    values$label_hjust <- as.integer(values$component != components[3])
    plot <- ggplot2::ggplot(values, ggplot2::aes(value, taxon, fill = component)) +
      ggplot2::geom_col(width = .65, colour = 'black', linewidth = .6,
        position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::scale_x_continuous(limits = c(0, limit), breaks = seq(0, limit, 10),
        expand = ggplot2::expansion(add = c(1.25, 3))) +
      ggplot2::geom_text(ggplot2::aes(x = label_x, y = label_y,
        label = sprintf('%.2f', value), colour = component, hjust = label_hjust),
        size = 5.1, fontface = 'bold', show.legend = FALSE) +
      ggplot2::geom_text(data = subset(values, component == components[3]),
        ggplot2::aes(x = separator_x, y = label_y), label = '|',
        inherit.aes = FALSE, colour = 'black', size = 4.8, fontface = 'bold') +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::scale_colour_manual(values = palette) +
      ggplot2::facet_wrap(ggplot2::vars(!!paste0('(', letters[i], ') ', forests[i]))) +
      ggplot2::labs(x = 'IVI points', y = NULL, fill = NULL) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(text = ggplot2::element_text(colour = 'black'),
        axis.text = ggplot2::element_text(size = 16, colour = 'black'),
        axis.text.y = ggplot2::element_text(face = 'italic'),
        axis.title.x = ggplot2::element_text(size = 18, face = 'bold'),
        axis.line.y = ggplot2::element_line(colour = 'black'),
        strip.background = ggplot2::element_rect(fill = 'grey90', colour = 'black'),
        strip.text = ggplot2::element_text(size = 17, face = 'bold', colour = 'black'),
        legend.position = 'inside', legend.position.inside = c(.98, .025),
        legend.justification = c(1, 0), legend.text = ggplot2::element_text(size = 16),
        legend.key.size = grid::unit(8, 'mm'),
        legend.background = ggplot2::element_rect(fill = 'white', colour = NA),
        plot.margin = ggplot2::margin(5.5, if (i == 1) 2.665 else 5.5,
          5.5, if (i == 2) 2.665 else 5.5))
    if (!nrow(z)) plot <- plot + ggplot2::annotate('text', x = limit/2, y = 1,
      label = 'Complete IVI unavailable')
    plot
  })
  print(patchwork::wrap_plots(panels, nrow = 1))
}

# ============================================================
# 10. Vs growth and ivi reporting
# ============================================================

# Workbook measurements are kept separate from the audited woody inventory.
vs_class_labels <- function(classes) {
  as.expression(lapply(classes, function(label) {
    if (startsWith(label, '>=')) {
      bquote('' >= .(substring(label, 3)))
    } else label
  }))
}

vs_draw_diameter <- function(counts, dbh_classes, forests, sample_sizes, heading, ylab) {
  forest_colours <- c(Ngangao = '#D4B58A', Mbololo = '#508A92')[forests]
  max_size <- 200
  ticks <- seq(0, max_size, 50)
  barplot(counts, beside = TRUE, names.arg = vs_class_labels(dbh_classes), col = forest_colours,
    xlab = 'Diameter class (cm)', ylab = ylab,
    ylim = c(-.035, 1.06) * max_size, yaxt = 'n', yaxs = 'i', cex.names = 16 / 12)
  title(heading, adj = 0, cex.main = 17.6 / 12)
  axis(2, at = ticks, las = 1)
  abline(h = par('usr')[3])
  legend('right', paste0(forests, ' (n = ', sample_sizes, ')'),
    pch = 22, pt.bg = forest_colours, pt.cex = 1.7, bty = 'n', cex = 16 / 12,
    inset = c((5 / 25.4) / par('pin')[1], 0))
}

vs_growth_and_ivi_reporting <- function(root, here, out, a, records, ivi, provisional = FALSE) {
  prefix <- if (provisional)'provisional_'else ''
  path <- function(n)rq_file(out, if (grepl('[.]md$', n))paste0(prefix, n)else n)
  put <- function(x, n)write.csv(x, path(paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  env <- new.env(parent = environment())
  source(file.path(here, '01_woody_composition.R'), local = env)
  ref <- env$floristic_reference(root)
  stages <- c('Listed larger stems', 'Recorded juveniles')
  forests <- c('Ngangao', 'Mbololo')
  key <- file.path(root, 'analysis/Datasets/derived/juveniles_record_taxon_counts.csv')
  if (!file.exists(key))env$run_floristic_composition(root, here)
  j <- read.csv(key, check.names = FALSE)
  stopifnot(setequal(j$record_id, records$record_id),
    all(j$forest == records$forest[match(j$record_id, records$record_id)]))
  x <- as.matrix(j[, setdiff(names(j), c('record_id', 'forest', 'status')), drop = FALSE])
  counts <- do.call(rbind, lapply(forests, function(f) {
        stem <- a[a$plot_id %in% records$record_id[records$forest == f] & a$woody_eligible,, drop = FALSE]
        larger <- table(stem$analysis_name)
        juvenile <- colSums(x[j$forest == f,, drop = FALSE], na.rm = TRUE)
        rbind(data.frame(forest = f, stage = stages[1], taxon = names(larger), count = as.numeric(larger)),
          data.frame(forest = f, stage = stages[2], taxon = names(juvenile), count = as.numeric(juvenile)))
      }))
  counts <- counts[counts$count>0,]
  labels <- c(T = 'Tree', S = 'Shrub', L = 'Liana')
  classes <- c(unname(labels), 'Multiple woody forms', 'Unresolved growth form')
  classify <- function(t) {
    h <- unique(ref$habit[which(ref$taxon == tolower(trimws(t)))])
    h <- h[!is.na(h) & nzchar(h)]
    if (!length(h) || !all(grepl('^[TSL](/[TSL])*$', h)))return(classes[5])
    h <- unique(unlist(strsplit(h, '/', fixed = TRUE)))
    if (length(h) == 1)unname(labels[h])else classes[4]
  }
  counts$growth_form <- vapply(counts$taxon, classify, character(1))
  put(counts, 'growth_form_taxon_count_audit')
  growth <- do.call(rbind, lapply(stages, function(stage)do.call(rbind, lapply(forests, function(f) {
            d <- counts[counts$stage == stage & counts$forest == f,]
            n <- vapply(classes, function(g)sum(d$count[d$growth_form == g]), numeric(1))
            surveyed <- if (stage == stages[1])sum(records$forest == f)else
            sum(j$forest == f & rowSums(!is.na(x))>0)
            data.frame(forest = f, stage, growth_form = classes, recorded_stems = unname(n),
              denominator_stems = sum(n), eligible_records = surveyed,
              share_percent = if (sum(n))100 * unname(n) / sum(n)else NA_real_)
          }))))
  put(growth, 'growth_form_proportions')
  primary <- ivi[ivi$forest %in% forests,]
  leading <- vs_table_s3(subset(primary, !is.na(rank) & rank <= 10))
  put(leading, 'Table_3_leading_taxa_by_forest')
  put(vs_table_s3(primary), 'Table_S5_Importance_Value_Index_among_listed_larger_stems')
  # The field protocol says above 5 DBH; the worksheets also list smaller stems.
  strict_stems <- a[which(is.finite(a$dbh) & a$dbh>5),]
  strict_ivi <- do.call(rbind, lapply(forests, function(f)
      vs_ivi(strict_stems, records$record_id[records$forest == f], f, 'DBH >5 cm sensitivity')))
  if (nrow(strict_ivi)) {
    compare <- merge(primary[, c('forest', 'taxon', 'listed_stems', 'rank', 'IVI_points')],
      strict_ivi[, c('forest', 'taxon', 'listed_stems', 'rank', 'IVI_points')],
      by = c('forest', 'taxon'), all.x = TRUE, suffixes = c('_listed', '_dbh_gt5'))
    compare$listed_stems_dbh_gt5[is.na(compare$listed_stems_dbh_gt5)] <- 0L
    compare$taxon <- sub('^Ocotea usambarensis$', 'Kuloa usambarensis',
      compare$taxon, ignore.case = TRUE)
    put(compare, 'IVI_protocol_DBH_threshold_sensitivity')
  }
  focus <- expand.grid(forest = forests, taxon = c('tabernaemontana stapfiana', 'albizia gummifera'),
    stringsAsFactors = FALSE)
  idx <- match(paste(focus$forest, focus$taxon), paste(primary$forest, tolower(primary$taxon)))
  focus$rank <- primary$rank[idx]
  focus$IVI_points <- primary$IVI_points[idx]
  focus$eligible_taxa_ranked <- vapply(focus$forest, function(f)sum(primary$forest == f & is.finite(primary$rank)), integer(1))
  put(focus, 'IVI_selected_taxa_by_forest')
  draw <- function() {
    par(mar = c(6.2, 5.1, 2.4, .65), oma = c(0, 0, 0, 0),
      mgp = c(2.85, .65, 0), cex = 1, cex.axis = 16 / 12, cex.lab = 18 / 12, font.lab = 2)
    pal <- c('#82A67B', '#DE9A54', '#999999', '#F4F4F4', '#CCCCCC')
    order_by_forest <- c(2, 4, 1, 3)
    bars <- matrix(growth$share_percent, nrow = length(classes))[, order_by_forest]
    n <- growth$denominator_stems[seq(1, nrow(growth), length(classes))][order_by_forest]
    mid <- barplot(bars, col = pal, border = 'black', space = c(.35, .65, 1.5, .65), ylim = c(-5, 115),
      xlim = c(0, 7.65), axisnames = FALSE, yaxt = 'n', yaxs = 'i',
      ylab = 'Share of eligible recorded stems (%)')
    title('Growth-form abundance', adj = 0, cex.main = 17.6 / 12)
    axis(2, at = seq(0, 100, 20), las = 1)
    axis(1, at = rowMeans(matrix(mid, ncol = 2, byrow = TRUE)), labels = rev(forests), line = 2.4,
      tick = FALSE, cex.axis = 16 / 12, font.axis = 2)
    abline(h = par('usr')[3])
    divider <- mean(mid[2:3])
    segments(divider, 0, divider, 111, lty = 3, col = 'grey55')
    text(mid, 109, as.expression(lapply(n, function(count) bquote(italic(n) == .(count)))),
      cex = 16 / 12)
    seedling_x <- c( - .12, - .12, - .52, - .93, - 1, - .64, - .26, 0, .26, .64, 1, .93, .52, .12, .12)
    seedling_y <- c( - 1, .02, .07, .38, .78, .73, .44, .16, .44, .73, .78, .38, .07, .02, - 1)
    for (i in seq_along(mid)) {
      x <- grconvertX(mid[i], from = 'user', to = 'inches')
      larger <- i %% 2L == 1L
      icon_scale <- .65
      icon_y <- icon_scale * (if (larger).15 else .105 * seedling_y)
      y <- grconvertY(par('usr')[3], from = 'user', to = 'inches') - 1 / 25.4 - max(icon_y)
      if (larger) {
        crown_x <- c( - .95, - .8, - .88, - .58, - .47, - .16, .05, .32, .48, .77, .74, 1, .86, .91, .55, .31, .05, - .19, - .47, - .68)
        crown_y <- c(.12, .35, .53, .59, .8, .77, 1, .93, .75, .72, .51, .36, .22, .02, - .06, .04, - .14, - .06, - .17, - .03)
        polygon(grconvertX(x + icon_scale * .15 * crown_x, 'inches', 'user'), grconvertY(y + icon_scale * .15 * crown_y, 'inches', 'user'),
          col = '#557747', border = NA, xpd = NA)
        polygon(grconvertX(x + icon_scale * ( - .035 + .105 * crown_x), 'inches', 'user'),
          grconvertY(y + icon_scale * (.035 + .105 * crown_y), 'inches', 'user'), col = '#789455', border = NA, xpd = NA)
        polygon(grconvertX(x + icon_scale * .15 * c( - .19, - .08, - .1, - .5, - .08, .04, .48, .15, .16, .25), 'inches', 'user'),
          grconvertY(y + icon_scale * .15 * c( - 1, - .58, - .14, .34, .08, .48, .44, .04, - .58, - 1), 'inches', 'user'),
          col = '#70513C', border = NA, xpd = NA)
      } else {
        polygon(grconvertX(x + icon_scale * .105 * seedling_x, 'inches', 'user'),
          grconvertY(y + icon_y, 'inches', 'user'), col = '#6F984C', border = NA, xpd = NA)
      }
      top <- cumsum(bars[, i])
      for (j in which(bars[, i]>0)) {
        share <- bars[j, i]
        label <- sprintf('%.1f%%', share)
        centre <- top[j] - share / 2
        if (share >= 5)text(mid[i], centre, label, col = 'black', font = 2, cex = 14 / 12)else {
          label_y <- centre + if (j == 2L) - 2 else 2
          segments(mid[i] + .5, centre, mid[i] + .59, label_y, lwd = 1.7)
          text(mid[i] + .59, label_y, label, pos = 4, offset = .1, cex = 14 / 12, font = 2)
        }
      }
    }
    legend_labels <- classes[1:4]
    key_width <- .075
    text_gap <- .03
    item_gap <- 2 / 25.4
    item_width <- key_width + text_gap + strwidth(legend_labels, units = 'inches', cex = 16 / 12)
    legend_centre <- mean(grconvertX(par('usr')[1:2], 'user', 'inches'))
    legend_left <- legend_centre - (sum(item_width) + 3 * item_gap) / 2
    key_left <- legend_left + c(0, cumsum(head(item_width, - 1) + item_gap))
    legend_y <- .11
    rect(grconvertX(key_left, 'inches', 'user'), grconvertY(legend_y - key_width / 2, 'inches', 'user'),
      grconvertX(key_left + key_width, 'inches', 'user'), grconvertY(legend_y + key_width / 2, 'inches', 'user'),
      col = pal[1:4], border = 'black', xpd = NA)
    text(grconvertX(key_left + key_width + text_gap, 'inches', 'user'), grconvertY(legend_y, 'inches', 'user'),
      legend_labels, adj = c(0, .5), cex = 16 / 12, xpd = NA)
  }
  for (ext in c('pdf', 'png')) {
    file <- path(paste0('Figure_S6_Growth_form_composition_by_fragment_and_size_class.', ext))
    if (ext == 'pdf')pdf(file, width = 12 / 2.54, height = 8.4 / 2.54,
      pointsize = 6, compress = dir.exists(tempdir()))
    else png(file, width = 12, height = 8.4, units = 'cm', res = 300, pointsize = 6)
    tryCatch(draw(), finally = dev.off())
  }

  invisible(growth)
}

# ============================================================
# 11. Vs provisional
# ============================================================

vs_provisional <- function(root, here, out, absent) {
  env <- new.env(parent = environment())
  source(file.path(here, '01_woody_composition.R'), local = env)
  input <- env$floristic_provisional(root, out, export = FALSE)
  rec <- input$records
  a <- input$source_rows[input$source_rows$stage == 'Listed larger stems',]
  a$included <- rec$included[match(a$record_id, rec$record_id)]
  number <- function(x) {
    x <- sub(',', '.', trimws(x), fixed = TRUE)
    suppressWarnings(as.numeric(ifelse(grepl('^[+]?[0-9]+([.][0-9]+)?$', x), x, NA_character_)))
  }
  a$dbh <- number(a$raw_dbh);
  a$height_m <- number(a$raw_height)
  a$dbh_decimal_comma <- !is.na(a$raw_dbh) & grepl('^[0-9]+,[0-9]+$', trimws(a$raw_dbh))
  a$height_decimal_comma <- !is.na(a$raw_height) & grepl('^[0-9]+,[0-9]+$', trimws(a$raw_height))
  a$basal_area_m2 <- vs_basal_area(a$dbh)
  a$plot_id <- a$record_id;
  a$analysis_name <- a$recorded_label
  habit_path <- 'analysis/Datasets/field/Data.Taita.xlsx'
  h <- readxl::read_excel(file.path(root, habit_path), sheet = 'Species.list', col_names = FALSE,
    col_types = 'text', .name_repair = 'minimal')
  norm <- function(x)tolower(trimws(gsub('[[:space:]\u00a0]+', ' ', x, perl = TRUE)))
  habits <- split(toupper(trimws(h[[4]])), norm(h[[3]]))
  a$source_habit <- vapply(a$recorded_label, function(x) {
      v <- unique(habits[[x]]);
    paste(v[!is.na(v) & nzchar(v)], collapse = '; ')
    }, character(1))
  a$woody_eligible <- grepl('^[TSL](/[TSL])*(; [TSL](/[TSL])*)*$', a$source_habit)
  a$habit_source <- habit_path
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  put(rec, 'worksheet_selection');
  put(a, 'stem_measurements')
  units <- input$evidence[input$evidence$stage == 'Listed larger stems',
    c('record_id', 'count_or_dbh_header', 'height_header')]
  put(units, 'worksheet_measurement_units')
  selected <- rec[rec$included,]
  med <- function(v) unname(vs_size_summary(v)['median'])
  records <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
        r <- selected[i,];
        z <- a[a$record_id == r$record_id & a$woody_eligible,]
        dbh <- z$dbh[is.finite(z$dbh) & z$dbh>0]
        d <- data.frame(record_id = r$record_id, forest = r$forest, listed_stems = nrow(z),
          dbh_available = length(dbh), height_available = sum(is.finite(z$height_m) & z$height_m>0),
          median_dbh_cm = med(z$dbh), median_height_m = med(z$height_m),
          quadratic_mean_dbh_cm = if (length(dbh))sqrt(mean(dbh^2))else NA_real_,
          observed_basal_area_m2 = if (length(dbh))sum(z$basal_area_m2, na.rm = TRUE)else NA_real_,
          basal_area_complete = length(dbh) == nrow(z),
          surveyed_area_ha = NA_real_, stems_per_ha = NA_real_, basal_area_m2_per_ha = NA_real_)
        h <- readxl::read_excel(file.path(root, r$workbook), sheet = r$sheet, range = 'A1:H10',
          col_names = FALSE, col_types = 'text', .name_repair = 'minimal')
        cells <- as.character(unlist(h))
        for (nm in c('canopy', 'shrub', 'herb', 'litter', 'bare')) {
          pattern <- switch (nm, canopy = 'upper canopy cover', shrub = 'cover of shrub layer',
            herb = 'cover of herbaceous layer', litter = 'cover of litter', bare = 'cover of bare')
          hit <- unique(cells[which(grepl(pattern, cells, ignore.case = TRUE))])
          v <- NA_real_
          if (length(hit) == 1L) {
            token <- regmatches(hit, regexpr(':[[:space:]]*[0-9]+([.,][0-9]+)?[[:space:]]*%', hit))
            if (length(token))v <- number(trimws(gsub('[:%]', '', token)))
          }
          d[[paste0(nm, '_percent')]] <- if (is.finite(v) && v <= 100)v else NA_real_
        }
        d
      }))
  put(records, 'structure_by_record')
  metrics <- c('listed_stems', 'median_dbh_cm', 'median_height_m', 'quadratic_mean_dbh_cm',
    'observed_basal_area_m2', 'canopy_percent', 'shrub_percent', 'herb_percent', 'litter_percent', 'bare_percent')
  summary <- do.call(rbind, lapply(split(records, records$forest), function(z)
      do.call(rbind, lapply(metrics, function(m) {
            v <- z[[m]][is.finite(z[[m]])]
            q <- if (length(v))unname(quantile(v, c(.25, .5, .75)))else rep(NA_real_, 3)
            data.frame(forest = z$forest[1], metric = m, survey_records = nrow(z), available_records = length(v),
              median = q[2], q1 = q[1], q3 = q[3])
          }))))
  summary$iqr <- summary$q3 - summary$q1
  units <- setNames(c('stems per record', 'cm', 'm', 'cm', 'm2 per record', rep('%', 5)), metrics)
  summary$unit <- unname(units[summary$metric])
  put(summary, 'forest_structure_summary')
  put(vs_table_2a(summary), 'Table_S4_Vegetation_structure_and_cover')
  contrast <- merge(summary[summary$forest == 'Mbololo', c('metric', 'unit', 'available_records', 'median')],
    summary[summary$forest == 'Ngangao', c('metric', 'unit', 'available_records', 'median')],
    by = c('metric', 'unit'), suffixes = c('_Mbololo', '_Ngangao'))
  contrast$median_difference <- with(contrast, median_Mbololo - median_Ngangao)
  put(contrast, 'structure_forest_contrasts')
  ivi <- do.call(rbind, lapply(split(selected, selected$forest), function(z)
      vs_ivi(a, z$record_id, z$forest[1], 'Worksheet labels')))
  ivi$status <- ifelse(is.na(ivi$rank), ivi$status,
    'IVI by worksheet label; taxonomy and stem eligibility await audit')
  put(ivi, 'worksheet_label_ivi')
  vs_growth_and_ivi_reporting(root, here, out, a, selected, ivi, provisional = TRUE)
  bins <- c(0, 5, 10, 20, 30, 50, Inf)
  classes <- c('<5', '5-<10', '10-<20', '20-<30', '30-<50', '>=50')
  dbh_counts <- do.call(rbind, lapply(c('Ngangao', 'Mbololo'), function(f) {
        z <- a[a$included & a$woody_eligible & a$forest == f,]
        cl <- cut(z$dbh, bins, labels = classes, right = FALSE)
        do.call(rbind, lapply(classes, function(k) {
              ix <- which(cl == k & is.finite(z$basal_area_m2))
              data.frame(forest = f, dbh_class = k, stems = length(ix), observed_basal_area_m2 = sum(z$basal_area_m2[ix]),
                valid_dbh = sum(is.finite(z$basal_area_m2)), missing_dbh = sum(!is.finite(z$basal_area_m2)))
            }))
      }))
  put(dbh_counts, 'diameter_class_distribution')
  save <- function(name, draw, w, h)for (ext in c('pdf', 'png')) {
    path <- rq_file(out, paste0(name, '.', ext))
    if (ext == 'pdf')pdf(path, width = w, height = h)else png(path, width = w * 200, height = h * 200, res = 200)
    tryCatch(draw(), finally = dev.off())
  }
  save('Figure_3_Vegetation_structure_and_cover', function() {
    vars <- c('listed_stems', 'median_dbh_cm', 'median_height_m', 'observed_basal_area_m2',
      'canopy_percent', 'shrub_percent', 'herb_percent', 'litter_percent')
    labels <- paste0('(', letters[seq_along(vars)], ') ', c('Stems per record',
      'Median DBH (cm)', 'Median height (m)', 'Basal area (m\u00b2/record)',
      'Upper canopy (%)', 'Shrub layer (%)', 'Herb layer (%)', 'Litter cover (%)'))
    forests <- c('Mbololo', 'Ngangao')
    values <- do.call(rbind, lapply(seq_along(vars), function(i) {
      data.frame(panel = factor(labels[i], levels = labels),
        forest = factor(records$forest, levels = forests), value = records[[vars[i]]])
    }))
    values <- values[is.finite(values$value) & !is.na(values$forest), ]
    groups <- split(values, list(values$panel, values$forest), drop = TRUE)
    stats <- do.call(rbind, lapply(groups, function(d) {
      q <- quantile(d$value, c(.25, .5, .75), names = FALSE)
      data.frame(panel = d$panel[1], forest = d$forest[1], q1 = q[1], median = q[2],
        q3 = q[3],
        sample_label = paste0('italic(n) == ', nrow(d)))
    }))
    violins <- do.call(rbind, Filter(function(d) length(unique(d$value)) > 1L, groups))
    if (is.null(violins)) violins <- values[FALSE, ]
    plot <- ggplot2::ggplot(values, ggplot2::aes(forest, value, fill = forest)) +
      ggplot2::geom_violin(data = violins, width = .72, scale = 'width', trim = TRUE,
        alpha = .48, colour = 'grey25') +
      ggplot2::geom_point(shape = 21, size = 3, stroke = .6, colour = 'black',
        position = ggplot2::position_jitter(width = .12, height = 0, seed = 42)) +
      ggplot2::geom_errorbar(data = stats, ggplot2::aes(forest, ymin = q1, ymax = q3),
        inherit.aes = FALSE, width = .16, linewidth = 1.1, colour = 'black') +
      ggplot2::geom_point(data = stats, ggplot2::aes(forest, median),
        inherit.aes = FALSE, size = 4.8, colour = 'black') +
      ggplot2::geom_text(data = stats, ggplot2::aes(forest, -Inf, label = sample_label),
        inherit.aes = FALSE, parse = TRUE, vjust = -.45, size = 5.5) +
      ggplot2::facet_wrap(~panel, ncol = 4, scales = 'free_y', drop = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.18, .18))) +
      ggh4x::facetted_pos_scales(y = list(
        panel == labels[5] ~ ggplot2::scale_y_continuous(
          limits = c(40, 100), breaks = seq(40, 100, 20),
          expand = ggplot2::expansion(add = c(6, 0))),
        panel == labels[6] ~ ggplot2::scale_y_continuous(
          limits = c(0, 60), breaks = seq(0, 60, 20), expand = c(0, 0))
      )) +
      ggplot2::scale_fill_manual(values = c(Ngangao = '#D4B58A', Mbololo = '#508A92')) +
      ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(), legend.position = 'bottom',
        legend.justification = 'center',
        legend.text = ggplot2::element_text(size = 17, colour = 'black'),
        legend.key.size = grid::unit(8, 'mm'),
        text = ggplot2::element_text(colour = 'black'),
        axis.text.y = ggplot2::element_text(colour = 'black', size = 15),
        strip.background = ggplot2::element_rect(fill = 'grey90'),
        strip.text = ggplot2::element_text(face = 'bold', size = 14, colour = 'black'),
        panel.grid.minor = ggplot2::element_blank())
    print(plot)
  }, 15, 7.5)
  save('Figure_4_Importance_Value_Index_by_fragment', function() vs_ranked_ivi_plot(ivi), 15, 7)
  save('Figure_S7_Stem_count_and_basal_area_by_DBH_class', function() {
    forests <- c('Mbololo', 'Ngangao')
    sample_sizes <- vapply(forests, function(f) sum(dbh_counts$stems[dbh_counts$forest == f]), numeric(1))
    d <- dbh_counts
    d$forest <- factor(d$forest, levels = forests)
    d$dbh_class <- factor(d$dbh_class, levels = classes)
    panels <- lapply(seq_len(2), function(i) {
      measure <- c('stems', 'observed_basal_area_m2')[i]
      maximum <- c(200, 20)[i]
      d$value <- d[[measure]]
      ggplot2::ggplot(d, ggplot2::aes(dbh_class, value, fill = forest)) +
        ggplot2::scale_x_discrete(labels = vs_class_labels) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = .8),
          width = .75, colour = 'black', linewidth = .45) +
        ggplot2::annotate('text', x = .58, y = maximum * 1.055,
          label = c('(a) Stem counts', '(b) Basal area')[i],
          hjust = 0, vjust = 1, fontface = 'bold', size = 6.2) +
        ggplot2::scale_y_continuous(breaks = seq(0, maximum, maximum/4),
          labels = function(x) sprintf('%.0f', x),
          expand = ggplot2::expansion(mult = c(.035, .075))) +
        ggplot2::coord_cartesian(ylim = c(0, maximum)) +
        ggplot2::scale_fill_manual(values = c(Ngangao = '#D4B58A', Mbololo = '#508A92'),
          labels = as.expression(Map(function(forest, count)
            bquote(.(forest) ~ '(' * italic(n) == .(count) * ')'), forests, sample_sizes))) +
        ggplot2::labs(x = 'Diameter class (cm)',
          y = if (i == 1) 'Listed stems with DBH' else expression(bold('Observed basal area (' * m^2 * ')')),
          fill = NULL) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(text = ggplot2::element_text(colour = 'black'),
          axis.title = ggplot2::element_text(face = 'bold', size = 18, colour = 'black'),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
          axis.text = ggplot2::element_text(size = 16, colour = 'black'),
          legend.text = ggplot2::element_text(size = 16, colour = 'black'),
          legend.key.size = grid::unit(8, 'mm'),
          legend.position = 'inside', legend.position.inside = c(if (i == 1) .97 else .04, .7),
          legend.justification = c(if (i == 1) 1 else 0, .5),
          legend.background = ggplot2::element_rect(fill = 'white', colour = NA))
    })
    print(patchwork::wrap_plots(panels, nrow = 1))
  }, 15, 7)
  status <- c('Provisional Question 2 completed. Audited results remain unavailable.',
    'Question 1 worksheet exclusions were used. Measurements and labels are from numbered larger-stem rows.',
    'Exact local Species.list matches with exclusively tree, shrub or liana habit codes were retained. No synonym decisions were inferred.',
    'DBH and height units were transcribed from each stem table header; see the worksheet measurement units audit. Measurement procedure and the accuracy of values cannot be reconstructed from those labels.',
    'Single decimal commas were parsed as decimal points and flagged in the stem audit. Only explicit percentage cover values were used.',
    'Verified surveyed area and independent resampling units are unavailable; per-hectare estimates and uncertainty intervals are withheld.',
    'Forest contrasts are Mbololo minus Ngangao in original units. Observed basal area may be partial.',
    'IVI = relative density (listed stem share) + relative frequency (share of summed taxon-record occurrences) + relative dominance (basal area share). Ordinary occurrence frequency uses survey records as its denominator.',
    'References: Wekesa et al. (2019), https://doi.org/10.1016/j.foreco.2019.03.011; Curtis and McIntosh (1951), https://doi.org/10.2307/1931725; Taubert et al. (2015), https://doi.org/10.1073/pnas.1513417112; Lutz et al. (2018), https://doi.org/10.1111/geb.12747.',
    'IVI ranks are withheld for forests with any missing or invalid DBH; worksheet labels have not been taxonomically reconciled.',
    'Missing audited inputs:', absent)

  message(paste(status, collapse = '\n'))
  invisible(list(records = records, summary = summary, ivi = ivi, output = out, provisional = TRUE))
}

# ============================================================
# 12. Run vegetation structure
# ============================================================

run_vegetation_structure <- function(root, here = file.path(root, 'analysis/scripts'), refresh = FALSE) {
  root <- vs_project_root(root);
  here <- file.path(root, 'analysis/scripts');
  vs_temp_ready()
  local_lib <- file.path(here, '.r-library')
  if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
  needed <- c('readxl', 'dplyr', 'jsonlite')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))stop('Missing packages: ', paste(missing, collapse = ', '),
    '. Install explicitly; this script never installs packages.')
  # Retained audited files live in the existing relevant analysis folders.
  source_paths <- rq_input(root, c('analysis/population_distribution/record_reconciliation.csv',
      'analysis/vegetation_structure/individual_stem_measurement_audit.csv',
      'analysis/floristic_composition/taxonomic_decisions.csv'))
  absent <- source_paths[!file.exists(source_paths)]
  out <- rq_dirs(root, 2)
  if (length(absent))return(vs_provisional(root, here, out, absent))
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  q <- read.csv(source_paths[1], check.names = FALSE)
  a <- read.csv(source_paths[2], check.names = FALSE)
  tax <- read.csv(source_paths[3], check.names = FALSE)
  stopifnot(!anyDuplicated(q$plot_id), !anyNA(match(a$plot_id, q$plot_id)),
    all(c('dbh', 'height_m', 'woody_eligible', 'included_primary') %in% names(a)))
  a$record_id <- a$plot_id
  a$source_stem_id <- paste(a$plot_id, a$sheet, a$source_row, sep = '::')
  if (anyDuplicated(a$source_stem_id))stop('Duplicate source-row identifier')
  a$dbh_status <- ifelse(is.na(a$dbh), 'Missing', ifelse(!is.finite(a$dbh),
      'Nonfinite', ifelse(a$dbh <= 0, 'Nonpositive', 'Positive recorded measurement')))
  a$height_status <- ifelse(is.na(a$height_m), 'Missing', ifelse(!is.finite(a$height_m),
      'Nonfinite', ifelse(a$height_m <= 0, 'Nonpositive', 'Positive recorded measurement')))
  a$basal_area_m2 <- vs_basal_area(a$dbh)
  a$measurement_units_verified <- FALSE
  a$dbh_class <- as.character(cut(a$dbh, c(0, 5, 10, 20, 30, 50, Inf), right = FALSE,
      labels = c('<5', '5-<10', '10-<20', '20-<30', '30-<50', '>=50')))
  a$dbh_class[a$dbh_status != 'Positive recorded measurement'] <- NA_character_
  a$height_class <- as.character(cut(a$height_m, c(0, 2, 5, 10, 20, Inf), right = FALSE,
      labels = c('<2', '2-<5', '5-<10', '10-<20', '>=20')))
  a$height_class[a$height_status != 'Positive recorded measurement'] <- NA_character_
  put(a, 'individual_stem_measurement_audit')
  put(q, 'record_eligibility');
  put(tax, 'taxonomic_decisions')
  write.csv(q, rq_input(root, 'record_eligibility.csv'), row.names = FALSE, na = 'NA')
  rows <- lapply(seq_len(nrow(q)), function(i) {
      z <- a[a$plot_id == q$plot_id[i] & a$woody_eligible,, drop = FALSE]
      ds <- vs_size_summary(z$dbh);
      hs <- vs_size_summary(z$height_m)
      good <- is.finite(z$dbh) & z$dbh>0
      d <- data.frame(record_id = q$plot_id[i], forest = q$forest[i], cohort = q$cohort[i],
        primary_eligible = q$primary_eligible[i], cultivated = q$cultivated[i],
        source_rows = sum(a$plot_id == q$plot_id[i]), listed_woody_stems = nrow(z),
        unresolved_or_excluded_woody_rows = sum(a$plot_id == q$plot_id[i] & !a$woody_eligible),
        dbh_valid_n = unname(ds['n']), dbh_missing = sum(is.na(z$dbh)),
        dbh_nonfinite = sum(!is.na(z$dbh) & !is.finite(z$dbh)),
        dbh_nonpositive = sum(is.finite(z$dbh) & z$dbh <= 0),
        height_valid_n = unname(hs['n']), height_missing = sum(is.na(z$height_m)),
        height_nonfinite = sum(!is.na(z$height_m) & !is.finite(z$height_m)),
        height_nonpositive = sum(is.finite(z$height_m) & z$height_m <= 0),
        quadratic_mean_dbh_cm = if (any(good))sqrt(mean(z$dbh[good]^2))else NA_real_,
        summed_basal_area_m2 = if (any(good))sum(z$basal_area_m2[good])else if (!nrow(z))0 else NA_real_,
        basal_area_measurements_n = sum(good), basal_area_complete = sum(good) == nrow(z))
      for (n in c('median', 'q1', 'q3', 'iqr', 'minimum', 'maximum')) {
        d[[paste0('dbh_', n, '_cm')]] <- unname(ds[n])
        d[[paste0('height_', n, '_m')]] <- unname(hs[n])
      }
      d
    })
  records <- dplyr::bind_rows(rows)
  put(records, 'vegetation_structure_by_record')
  cohorts <- list(primary = q$plot_id[q$primary_eligible],
    broader_forest_context = q$plot_id[!q$cultivated], cultivated = q$plot_id[q$cultivated])
  ivi_rows <- list();
  ii <- 0L
  for (cohort in c('primary', 'broader_forest_context'))for (forest in c('Ngangao', 'Mbololo')) {
    ids <- intersect(cohorts[[cohort]], q$plot_id[q$forest == forest])
    ii <- ii + 1L;
    ivi_rows[[ii]] <- vs_ivi(a, ids, forest, cohort)
  }
  ivi <- dplyr::bind_rows(ivi_rows)
  put(ivi, 'importance_value_index')
  top_ivi <- ivi[ivi$cohort == 'primary' & !is.na(ivi$rank) & ivi$rank <= 10,]
  put(top_ivi, 'importance_value_index_top_taxa')
  primary_records <- data.frame(record_id = q$plot_id[q$primary_eligible], forest = q$forest[q$primary_eligible])
  vs_growth_and_ivi_reporting(root, here, out, a, primary_records, ivi[ivi$cohort == 'primary',])
  ivi_check <- dplyr::bind_rows(lapply(split(ivi, interaction(ivi$cohort, ivi$forest, drop = TRUE)), function(z)
      data.frame(cohort = z$cohort[1], forest = z$forest[1], survey_records = z$survey_records[1],
        taxon_labels = nrow(z), listed_stems = z$forest_listed_stem_total[1],
        taxon_record_occurrences = z$summed_taxon_record_occurrences[1],
        observed_basal_area_m2 = z$forest_observed_basal_area_m2[1],
        missing_or_invalid_dbh = sum(z$dbh_missing_or_invalid),
        relative_abundance_sum = sum(z$relative_stem_abundance_percent),
        relative_frequency_sum = sum(z$relative_record_frequency_percent),
        relative_dominance_sum = sum(z$relative_dominance_percent),
        IVI_points_sum = sum(z$IVI_points), IVI_percent_sum = sum(z$IVI_percent))))
  put(ivi_check, 'importance_value_index_denominators')
  # Forest summaries give each record one observation. Individual stem counts
  # and size measurements are never used as independent forest replicates.
  metrics <- c('listed_woody_stems', 'dbh_median_cm', 'dbh_iqr_cm', 'quadratic_mean_dbh_cm',
    'summed_basal_area_m2', 'height_median_m', 'height_iqr_m', 'height_maximum_m')
  summaries <- list();
  k <- 0L
  for (cohort in names(cohorts))for (forest in c('Ngangao', 'Mbololo')) {
    z <- records[records$record_id %in% cohorts[[cohort]] & records$forest == forest,]
    if (!nrow(z))next
    for (metric in metrics) {
      v <- z[[metric]];
    valid <- v[is.finite(v)]
      qs <- if (length(valid))unname(quantile(valid, c(.25, .5, .75), type = 7))else rep(NA_real_, 3)
      k <- k + 1L;
      summaries[[k]] <- data.frame(cohort, forest, metric,
        survey_records = nrow(z), available_records = length(valid), missing_records = sum(!is.finite(v)),
        median = qs[2], q1 = qs[1], q3 = qs[3], iqr = qs[3] - qs[1],
        minimum = if (length(valid))min(valid)else NA_real_,
        maximum = if (length(valid))max(valid)else NA_real_,
        observation_unit = 'Source survey record; forest representation and independence unverified')
    }
  }
  summary <- dplyr::bind_rows(summaries);
  put(summary, 'forest_record_summary')
  # Report all common size classes, including zero counts, with the correct
  # record-specific denominator separately for diameter and height.
  bins <- list(dbh = c('<5', '5-<10', '10-<20', '20-<30', '30-<50', '>=50'),
    height = c('<2', '2-<5', '5-<10', '10-<20', '>=20'))
  class_rows <- list();
  k <- 0L
  for (i in seq_len(nrow(q)))for (dimension in names(bins)) {
    z <- a[a$plot_id == q$plot_id[i] & a$woody_eligible,, drop = FALSE]
    classes <- z[[paste0(dimension, '_class')]];
  den <- sum(!is.na(classes))
    counts <- as.integer(table(factor(classes, levels = bins[[dimension]])))
    k <- k + 1L;
    class_rows[[k]] <- data.frame(record_id = q$plot_id[i], forest = q$forest[i],
      primary_eligible = q$primary_eligible[i], dimension, class = bins[[dimension]],
      stems = counts, measurement_denominator = den,
      proportion = if (den)counts / den else NA_real_,
      missing_or_invalid = nrow(z) - den)
  }
  classes <- dplyr::bind_rows(class_rows);
  put(classes, 'size_classes_by_record')
  pooled <- list();
  k <- 0L
  for (forest in c('Ngangao', 'Mbololo'))for (dimension in names(bins)) {
    z <- classes[classes$primary_eligible & classes$forest == forest & classes$dimension == dimension,]
    ids <- unique(z$record_id);
    den <- sum(z$measurement_denominator[!duplicated(z$record_id)])
    counts <- vapply(bins[[dimension]], function(cl)sum(z$stems[z$class == cl]), numeric(1))
    k <- k + 1L;
    pooled[[k]] <- data.frame(forest, dimension, class = bins[[dimension]],
      survey_records = length(ids), pooled_stems = counts, measurement_denominator = den,
      proportion = if (den)counts / den else NA_real_,
      interpretation = 'Supplementary pooled stems; records contribute unequally')
  }
  pooled <- dplyr::bind_rows(pooled);
  put(pooled, 'supplementary_pooled_size_distributions')
  mapping <- c(canopy = 'Upper canopy cover', shrub_cover = 'Shrub layer cover',
    herb_cover = 'Herbaceous layer cover', litter = 'Litter cover', bare = 'Bare ground and rocks')
  cover <- dplyr::bind_rows(lapply(names(mapping), function(field) {
        v <- q[[field]]
        data.frame(record_id = q$plot_id, forest = q$forest, primary_eligible = q$primary_eligible,
          attribute = unname(mapping[field]), domain = if (field %in% c('litter', 'bare'))
          'Ground surface'else 'Vegetation layer', recorded_percent = v,
          valid = is.finite(v) & v >= 0 & v <= 100,
          status = ifelse(is.na(v), 'Missing', ifelse(!is.finite(v), 'Nonfinite',
              ifelse(v<0 | v>100, 'Outside 0 to 100 percent', 'Recorded cover'))))
      }))
  put(cover, 'vegetation_layers_and_ground_cover_by_record')
  cs <- list();
  k <- 0L
  for (forest in c('Ngangao', 'Mbololo'))for (attribute in unname(mapping)) {
    z <- cover[cover$primary_eligible & cover$forest == forest & cover$attribute == attribute,]
    v <- z$recorded_percent[z$valid];
    qs <- if (length(v))unname(quantile(v, c(.25, .5, .75)))else rep(NA_real_, 3)
    k <- k + 1L;
    cs[[k]] <- data.frame(forest, attribute, domain = z$domain[1],
      survey_records = nrow(z), available_records = length(v), missing_or_invalid = sum(!z$valid),
      median_percent = qs[2], q1_percent = qs[1], q3_percent = qs[3], iqr_percent = qs[3] - qs[1],
      minimum_percent = if (length(v))min(v)else NA_real_, maximum_percent = if (length(v))max(v)else NA_real_)
  }
  put(dplyr::bind_rows(cs), 'forest_cover_summary')
  # Separate missing, invalid and available measurements, including records with
  # no usable height. No missing size is replaced with zero or estimated.
  put(records[, c('record_id', 'forest', 'primary_eligible', 'listed_woody_stems',
        'dbh_valid_n', 'dbh_missing', 'dbh_nonfinite', 'dbh_nonpositive',
        'height_valid_n', 'height_missing', 'height_nonfinite', 'height_nonpositive',
        'basal_area_complete')], 'measurement_completeness')
  cols <- c(Ngangao = '#D4B58A', Mbololo = '#508A92')
  panel <- function(value, forest, label, heading) {
    ok <- is.finite(value)
    f <- factor(forest, levels = names(cols))
    if (!any(ok)) {plot.new();
      title(main = heading);
      text(.5, .5, 'No usable measurements');
      return(invisible())}
    groups <- split(ifelse(ok, value, NA_real_), f, drop = FALSE)
    limits <- range(value[ok]);
    padding <- if (diff(limits)>0)diff(limits) * .06 else max(1, abs(limits[1])) * .06
    # Explicit limits retain every record, including extreme values beyond the
    # boxplot whiskers. outline=FALSE only prevents duplicated outlier symbols.
    boxplot(groups, col = adjustcolor(cols, alpha.f = .18),
      names = paste(names(cols), paste0('n=', vapply(groups, function(z)sum(is.finite(z)), integer(1))), sep = '\n'),
      ylim = limits + c( - padding, padding), outline = FALSE, ylab = label, xlab = '', main = heading, cex.main = .9)
    # Deterministic offsets show every record without random jitter.
    for (i in seq_along(cols)) {
      ix <- which(ok & as.integer(f) == i)
      if (!length(ix))next
      off <- if (length(ix)>1)seq( - .11, .11, length.out = length(ix))else 0
      points(i + off, value[ix], pch = 16, col = cols[i], cex = .7)
    }
  }
  main_metrics <- c('listed_woody_stems', 'dbh_median_cm', 'quadratic_mean_dbh_cm',
    'summed_basal_area_m2', 'height_median_m', 'height_iqr_m')
  main_labels <- c('Listed stems per record', 'Median DBH (cm)', 'Quadratic mean DBH (cm)',
    'Summed basal area per record (m2)', 'Median listed-stem height (m)', 'Within-record height IQR (m)')
  draw_structure <- function() {
    par(mfrow = c(2, 3), mar = c(4, 4.4, 3.3, 1))
    z <- records[records$primary_eligible,]
    for (i in seq_along(main_metrics))panel(z[[main_metrics[i]]], z$forest, main_labels[i],
      paste0(letters[i], ') ', c('Listed abundance', 'Diameter', 'Quadratic mean diameter',
          'Listed basal area', 'Height', 'Height variation')[i]))
  }
  draw_cover <- function() {
    par(mfrow = c(2, 3), mar = c(4, 4.4, 3.3, 1))
    for (i in seq_along(mapping)) {
      z <- cover[cover$primary_eligible & cover$attribute == mapping[i],]
      panel(ifelse(z$valid, z$recorded_percent, NA_real_), z$forest, 'Recorded cover (%)',
        paste0(letters[i], ') ', mapping[i]))
    }
    plot.new();
    text(.5, .6, 'Vegetation layers and ground surfaces\nare separate attributes.', cex = .9)
    text(.5, .35, 'Each point represents one source record.\nNo composite score is calculated.', cex = .85)
  }
  draw_pooled <- function() {
    par(mfrow = c(2, 2), mar = c(5, 4.4, 3.3, 1))
    for (dimension in names(bins))for (forest in names(cols)) {
      z <- pooled[pooled$dimension == dimension & pooled$forest == forest,]
      barplot(z$proportion, names.arg = vs_class_labels(z$class), col = cols[forest], las = 2,
        ylim = c(0, 1), ylab = 'Proportion of pooled valid measurements',
        main = paste(forest, if (dimension == 'dbh')'DBH (cm)'else 'height (m)'), cex.main = .9)
      mtext(paste('Measurement denominator =', z$measurement_denominator[1]), side = 3, line = .3, cex = .7)
    }
  }
  draw_ivi <- function()vs_ranked_ivi_plot(top_ivi)
  draws <- list(vegetation_structure_record_comparisons = draw_structure,
    vegetation_layers_and_ground_cover = draw_cover,
    supplementary_pooled_size_distributions = draw_pooled,
    importance_value_index_components = draw_ivi)
  figure_names <- c(vegetation_structure_record_comparisons = 'Figure_2_recorded_vegetation_structure',
    vegetation_layers_and_ground_cover = 'Figure_S2_vegetation_layers_and_ground_cover',
    supplementary_pooled_size_distributions = 'Figure_S3_pooled_diameter_and_height_distributions',
    importance_value_index_components = 'Figure_S4_importance_value_index_components')
  for (name in names(draws))for (ext in c('pdf', 'png')) {
    vs_temp_ready()
    if (ext == 'pdf')pdf(rq_file(out, paste0(figure_names[[name]], '.pdf')), width = 12, height = 8)
    else png(rq_file(out, paste0(figure_names[[name]], '.png')), width = 3000, height = 2000, res = 250)
    tryCatch(draws[[name]](), finally = dev.off())
  }
  notes <- c('# Vegetation structure methods', '',
    'How does vegetation structure differ between the surveyed neighbourhoods of the two fragments, and which woody taxa have the highest IVI in each fragment?',
    'This analysis reads the retained population-distribution record ledger, individual-stem measurement audit and floristic taxonomic decisions from their existing analysis folders. Structural metrics and IVI are recalculated from audited stem measurements; source workbooks are not re-extracted. No configuration or external reporting-helper file is required. Primary summaries use 23 forest-associated records: nine Ngangao and 14 Mbololo. Cultivated cases and identity conflicts remain separately audited. Reference-mapped woody eligibility and source-name decisions are preserved. Headers are never added as extra measured stems.',
    'Each listed source row is a stem observation, not a verified unique plant or genet. All extracted rows remain in individual_stem_measurement_audit.csv; structural metrics use eligible woody rows only. Source extraction includes associated-stem rows with a numeric DBH entry, so the audit cannot certify completeness of unextracted or unrecorded plants.',
    'Positive finite DBH and height are summarised independently. Type-7 quantiles provide median, Q1, Q3 and IQR. Quadratic mean DBH is sqrt(mean(DBH^2)). Individual basal area is pi*(DBH_cm/200)^2 in m2, then summed per record; diameters are never summed before squaring. The BA denominator and completeness flag identify partial measurement sums. No usable diameter in a nonempty inventory gives NA, while an empty eligible inventory has zero listed basal area.',
    'Missing, nonfinite and nonpositive measurements remain distinct. They are not imputed or set to zero. Reported units follow source labels and are not independently certified; unusual measurements are retained for review. Maximum listed-stem height is not labelled canopy height.',
    'DBH classes are <5, 5-<10, 10-<20, 20-<30, 30-<50 and >=50 cm; height classes are <2, 2-<5, 5-<10, 10-<20 and >=20 m. Each record retains a separate valid measurement denominator. Pooled distributions are supplementary because records contribute different numbers of stems.',
    'Forest summaries are medians, quartiles, IQR and ranges of record-level metrics, with available and missing record counts. Records receive equal weight. These summaries do not estimate the median size of all trees in either forest. Record points are plotted with deterministic horizontal offsets; all extremes are retained within explicit axis limits. Boxes use R boxplot Tukey hinges and whiskers extending to observations within 1.5 hinge spreads; tabulated quartiles use type 7. Labels show available record n. No inferential tests or confidence limits are shown.',
    'Upper canopy, shrub and herbaceous cover remain separate vegetation attributes. Litter and bare ground/rock cover are ground-surface attributes. Missing or out-of-range values are flagged and excluded from cover summaries; raw values remain in the audit. Cover values are never added into a complexity or disturbance score.',
    'Importance Value Index (IVI) is calculated separately within each forest and cohort for listed woody stems, using accepted analysis labels. Relative stem abundance is 100*n_i/sum(n_i); relative record frequency is 100*f_i/sum(f_i), where f_i is the number of source records listing taxon i; relative dominance is 100*BA_i/sum(BA_i), based on summed individual basal areas. The frequency denominator is summed taxon-record occurrences, not number of surveyed records. Ordinary record frequency (100*f_i/number of records) is also reported separately.',
    'IVI points sum the three percentage components (total 300 across all eligible taxa within a forest/cohort); IVI_percent divides points by three (total 100). This is a unitless equal-weight composite, not an area-based density estimate. Primary and broader forest-context alternatives remain separate; cultivated records, juvenile aggregate counts and hypothetical focal-header additions do not enter IVI. Missing or invalid DBH makes complete dominance and IVI undefined for the entire forest/cohort, while observed basal-area shares remain audit quantities. Ranks use minimum rank for ties; the figure includes all taxa at ranks up to 10.',
    'IVI describes structural prominence in selected source records, not forest-wide dominance, conservation priority, reproductive importance or disturbance tolerance. Unequal footprints, focal-targeted placement, possible focal omission, unverified measurement units and field identity can affect rankings. No inferential comparison of IVI ranks is performed.',
    'No per-hectare scaling, biomass, carbon, demographic trend, recovery, successional trajectory or disturbance effect is estimated. Unknown footprint, inclusion thresholds, observer effects, measurement units, overlap and representativeness limit regional inference.',
    '', 'Sources and scope',
    'Barandun et al. (2026), Journal of Ecology 114, e70232: sampling area and observer comparability. https://doi.org/10.1111/1365-2745.70232',
    'Matsuo et al. (2021), Journal of Ecology 109, 2871-2884: individual stem measurements and plot-level structural attributes. https://doi.org/10.1111/1365-2745.13680',
    'Coverdale and Davies (2023), Journal of Ecology 111, 1378-1395: complementary structural attributes rather than one universal index. https://doi.org/10.1111/1365-2745.14068',
    'These descriptive calculations are adaptations to the available records; the published sampling designs are not reproduced.',
    '', 'Run: Rscript --vanilla "analysis/scripts/02_vegetation_structure.R"',
    'Outputs: analysis/outputs/research question 2 (figures/ and tables/). No manuscript is edited.')
  writeLines(notes, rq_file(out, 'METHODS.md'))
  results <- c('# Research question 2: Vegetation structure and IVI', '',
    'Primary summaries describe survey records, not complete regional vegetation.',
    '', '## Available measurements')
  for (f in names(cols)) {
    z <- records[records$primary_eligible & records$forest == f,]
    results <- c(results, sprintf('%s: %d records, %d listed eligible woody stems; %d positive DBH measurements; %d positive heights and %d missing heights.',
        f, nrow(z), sum(z$listed_woody_stems), sum(z$dbh_valid_n), sum(z$height_valid_n), sum(z$height_missing)))
  }
  results <- c(results, '', '## Record-level forest summaries')
  for (i in which(summary$cohort == 'primary')) {
    z <- summary[i,]
    results <- c(results, sprintf('%s, %s: median %.3f; Q1-Q3 %.3f-%.3f; %d/%d available records.',
        z$forest, z$metric, z$median, z$q1, z$q3, z$available_records, z$survey_records))
  }
  paired_summary <- summary[summary$cohort == 'primary', c('forest', 'metric', 'median', 'available_records', 'survey_records')]
  contrast <- merge(paired_summary[paired_summary$forest == 'Mbololo',],
    paired_summary[paired_summary$forest == 'Ngangao',], by = 'metric', suffixes = c('_Mbololo', '_Ngangao'))
  contrast$median_difference <- with(contrast, median_Mbololo - median_Ngangao)
  put(contrast, 'primary_structure_forest_contrasts')
  results <- c(results, '', 'Descriptive differences in record medians (Mbololo minus Ngangao) are in primary_structure_forest_contrasts.csv.',
    '', '## Importance Value Index')
  for (forest in names(cols)) {
    z <- ivi[ivi$cohort == 'primary' & ivi$forest == forest & !is.na(ivi$rank),]
    z <- z[order(z$rank, z$taxon),]
    if (nrow(z)) {
      top <- head(z, 3)
      results <- c(results, paste0(forest, ': ', paste(sprintf('%s (%.2f IVI points)', top$taxon, top$IVI_points), collapse = '; '), '.'))
    } else results <- c(results, paste0(forest, ': complete IVI unavailable; see denominator audit.'))
  }
  results <- c(results, 'IVI rankings describe listed-stem structural prominence only; they do not establish conservation importance or whole-forest dominance.')
  writeLines(results, rq_file(out, 'RESULTS.md'))
  export_vegetation_reporting(out)
  inputs <- c(file.path(here, '02_vegetation_structure.R'), source_paths)
  put(data.frame(path = substring(inputs, nchar(root) + 2), md5 = unname(tools::md5sum(inputs))),
    'input_and_code_hashes')
  writeLines(capture.output(sessionInfo()), rq_file(out, 'sessionInfo.txt'))
  jsonlite::write_json(list(primary_records = sum(q$primary_eligible),
      observation_unit = 'Source survey record', dbh_unit = 'Recorded centimetres',
      height_unit = 'Recorded metres', basal_area_unit = 'm2 per record',
      area_verified = FALSE, measurement_units_verified = FALSE, confirmatory_tests = FALSE,
      importance_value_index = TRUE, IVI_components = c('Relative listed-stem abundance',
        'Relative record frequency', 'Relative individual-summed basal area'),
      forest_summary = 'Distribution of record-level metrics; records weighted equally'),
    rq_file(out, 'analysis_status.json'), auto_unbox = TRUE, pretty = TRUE)
  message('Vegetation structure complete: ', out)
  invisible(list(records = records, summary = summary, cover = cover, output = out))
}

# ============================================================
# 13. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep('^--file=', commandArgs(), value = TRUE)
  start <- if (length(arg))dirname(normalizePath(sub('^--file=', '', arg[1]), mustWork = TRUE))else getwd()
  root <- vs_project_root(start)
  run_vegetation_structure(root, file.path(root, 'analysis/scripts'))
  source(file.path(file.path(root, 'analysis'), 'scripts', 'combine_output_tables.R'), local = TRUE)
  combine_output_tables(file.path(file.path(root, 'analysis'), 'outputs', 'research question 2'), '02_vegetation_structure.R')
}
