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
# 1. Jao project root
# ============================================================

# Research question 4: For each focal species, what proportion of all survey records in each fragment contains both juveniles and larger stems, and how does this proportion differ between Ngangao and Mbololo?
# Descriptive juvenile abundance, record frequency and stage co-occurrence.

jao_project_root <- function(start = getwd()) {
  if (length(start) != 1L || is.na(start) || !nzchar(start))
  stop('Provide an existing project directory.', call. = FALSE)
  p <- normalizePath(start, winslash = '/', mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, 'analysis/scripts/04_size_class_overlap.R')))return(p)
    parent <- dirname(p)
    if (identical(parent, p))break
    p <- parent
  }
  stop('Cannot locate the project from ', start,
    '. Set the working directory to the project or a folder within it.', call. = FALSE)
}

# ============================================================
# 2. Jao category
# ============================================================

jao_category <- function(j, l) {
  out <- rep(NA_character_, length(j))
  ok <- is.finite(j) & is.finite(l)
  out[ok] <- ifelse(j[ok]>0, ifelse(l[ok]>0, 'Both stages', 'Juveniles only'),
    ifelse(l[ok]>0, 'Larger stems only', 'Neither stage listed'))
  out
}

# ============================================================
# 3. Jao summarise
# ============================================================

jao_summarise <- function(z) {
  x <- z$juvenile_count[is.finite(z$juvenile_count)]
  paired <- !is.na(z$cooccurrence)
  cats <- c('Both stages', 'Juveniles only', 'Larger stems only', 'Neither stage listed')
  counts <- vapply(cats, function(k)sum(z$cooccurrence == k, na.rm = TRUE), integer(1))
  qs <- if (length(x))unname(quantile(x, c(.25, .5, .75), type = 7))else rep(NA_real_, 3)
  denom <- counts[1] + counts[2]
  abundance <- data.frame(source_records = nrow(z), eligible_juvenile_records = length(x),
    missing_juvenile_records = sum(!is.finite(z$juvenile_count)),
    juvenile_total = if (length(x))sum(x)else NA_real_, median = qs[2], q1 = qs[1], q3 = qs[3],
    minimum = if (length(x))min(x)else NA_real_, maximum = if (length(x))max(x)else NA_real_,
    juvenile_positive_records = sum(x>0), record_frequency = if (length(x))mean(x>0)else NA_real_,
    record_frequency_percent = if (length(x))100 * mean(x>0)else NA_real_,
    paired_eligible_records = sum(paired), missing_either_stage_records = sum(!paired),
    paired_juvenile_positive_records = unname(denom),
    both_stages_records = unname(counts[1]),
    larger_presence_given_juvenile = if (denom>0)unname(counts[1] / denom)else NA_real_)
  categories <- data.frame(category = cats, records = unname(counts), paired_eligible_records = sum(paired),
    percent = if (sum(paired)>0)100 * unname(counts) / sum(paired)else NA_real_)
  list(abundance = abundance, categories = categories)
}

# ============================================================
# 4. Jao overlap report
# ============================================================

jao_overlap_report <- function(z, out, provisional = FALSE) {
  forests <- c('Ngangao', 'Mbololo')
  species <- c('coffea fadenii', 'afrocarpus usambarensis')

  category_labels <- c('Both', 'Juveniles only', 'Larger stems only', 'Neither')
  category_colours <- c('#297A3F', '#E9CE63', '#8BB657', '#E5E5E5')
  output <- function(name)rq_file(out, if (provisional && grepl('[.]md$', name))paste0('provisional_', name)else name)
  put <- function(x, name)write.csv(x, output(paste0(name, '.csv')), row.names = FALSE, na = 'NA')
  stopifnot(!anyDuplicated(z[c('record_id', 'species')]))
  z$joint_detection <- ifelse(is.na(z$cooccurrence), NA_integer_,
    as.integer(z$cooccurrence == 'Both stages'))
  groups <- expand.grid(forest = forests, species = species, stringsAsFactors = FALSE)
  summaries <- categories <- vector('list', nrow(groups))
  for (i in seq_len(nrow(groups))) {
    g <- groups[i,]
    d <- z[z$forest == g$forest & z$species == g$species,]
    n <- sum(!is.na(d$joint_detection))
    both <- sum(d$joint_detection, na.rm = TRUE)
    summaries[[i]] <- data.frame(g, selected_records = nrow(d), paired_eligible_records = n,
      missing_either_class_records = nrow(d) - n, both_size_classes_records = both,
      recorded_overlap_percent = if (n)100 * both / n else NA_real_)
    categories[[i]] <- cbind(g[rep(1, 4),], jao_summarise(d)$categories)
  }
  overlap <- do.call(rbind, summaries)
  categories <- do.call(rbind, categories)
  contrast <- merge(overlap[overlap$forest == forests[2], - 1],
    overlap[overlap$forest == forests[1], - 1], by = 'species',
    suffixes = c('_mbololo', '_ngangao'), sort = FALSE)
  contrast$difference_mbololo_minus_ngangao_pp <- with(contrast,
    recorded_overlap_percent_mbololo - recorded_overlap_percent_ngangao)
  put(overlap, 'recorded_size_class_overlap_summary')
  manuscript <- unique(categories[c('species', 'forest', 'paired_eligible_records')])
  category_names <- c('Both stages' = 'both', 'Juveniles only' = 'juveniles_only',
    'Larger stems only' = 'listed_larger_stems_only', 'Neither stage listed' = 'neither')
  for (category in names(category_names)) {
    values <- categories[categories$category == category,]
    index <- match(paste(manuscript$species, manuscript$forest), paste(values$species, values$forest))
    manuscript[[category_names[[category]]]] <- values$records[index]
  }
  put(manuscript, 'Table_4_recorded_size_class_overlap')
  put(contrast, 'size_class_overlap_forest_contrasts')
  put(categories, 'four_category_summary')
  put(z, 'joint_detection_record_audit')
  put(z[is.na(z$joint_detection),], 'missing_size_class_records')
  draw_categories <- function() {
    category_order <- c(1L, 3L, 2L, 4L)
    par(mfrow = c(1, 2), mar = c(3, 3.4, 2.4, .65), oma = c(3.2, 2.1, 0, 0),
      mgp = c(2.85, .65, 0), cex = 1, cex.axis = 17/12)
    forest_order <- rev(forests)
    for (i in seq_along(species)) {
      d <- categories[categories$species == species[i], ]
      values <- matrix(d$percent, nrow = 4)[category_order, match(forest_order, forests), drop = FALSE]
      mid <- barplot(values, col = category_colours[category_order], border = 'black',
        ylim = c(-5, 115), yaxt = 'n', yaxs = 'i', axisnames = FALSE)
      axis(2, at = seq(0, 100, 20), las = 1)
      axis(1, at = mid, labels = forest_order, tick = FALSE, font.axis = 2)
      abline(h = par('usr')[3])
      sp <- species[i]
      name <- paste0(toupper(substr(sp, 1, 1)), substring(sp, 2))
      title(bquote(bold(.(paste0('(', letters[i], ')'))) ~ italic(.(name))),
        adj = 0, cex.main = 18.6/12)
      n <- overlap[overlap$species == sp, ]
      n <- n$paired_eligible_records[match(forest_order, n$forest)]
      text(mid, 108, as.expression(lapply(n, function(count) bquote(italic(n) == .(count)))),
        cex = 17/12)
      for (j in seq_along(mid)) {
        positive <- which(values[, j] > 0)
        centres <- cumsum(values[, j]) - values[, j]/2
        text(mid[j], centres[positive], sprintf('%.1f%%', values[positive, j]),
          font = 2, cex = 16/12)
      }
    }
    mtext('Paired survey records (%)', side = 2, outer = TRUE, line = .5, cex = 19/12, font = 2)
    par(fig = c(0, 1, 0, .105), new = TRUE, mar = rep(0, 4), oma = rep(0, 4))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
    keys <- c(.18, .35, .57, .75)
    for (position in seq_along(keys)) {
      i <- category_order[position]
      x <- keys[position]
      rect(x, .34, x + .012, .61, col = category_colours[i], border = 'black')
      if (i %in% c(1, 3)) {
        tree <- x - if (i == 1L) .034 else .015
        rect(tree - .0017, .19, tree + .0017, .51, col = '#70513C', border = NA)
        polygon(tree + 1.3 * c(-.008, -.006, -.007, -.003, .001, .006, .008, .005, -.004),
          .5 + 1.3 * (c(.48, .59, .69, .76, .73, .76, .59, .48, .44) - .5),
          col = '#557747', border = NA)
      }
      if (i %in% c(1, 2)) {
        seed <- x - .013
        polygon(seed + .0091*c(-.12, -.12, -.52, -1, -.64, 0, .64, 1, .52, .12, .12),
          .48 + .247*c(-1, 0, .1, .8, .72, .16, .72, .8, .1, 0, -1),
          col = '#6F984C', border = NA)
      }
      if (i == 4L) text(x - .012, .47, '-', cex = 1.2)
      text(x + .017, .47, category_labels[i], adj = c(0, .5), cex = 17/12)
    }
  }
  for (ext in c('pdf', 'png')) {
    path <- output(paste0('Figure_6_Co_occurrence_of_juveniles_and_listed_larger_stems.', ext))
    if (ext == 'pdf') cairo_pdf(path, width = 15, height = 7, family = 'Arial')
    else png(path, width = 3300, height = 1540, res = 220, type = 'cairo', family = 'Arial')
    tryCatch(draw_categories(), finally = dev.off())
  }
  methods <- c('# Research question 4. Recorded size-class overlap', '',
    if (provisional)'This run uses the provisional Question 3 worksheet ledger and its primary record selection.'else
    'This run uses the primary audited forest-associated records.',
    'Each record contributes one binary joint-detection value per focal species. Both counts must be available. Positive counts of both classes give 1; every other paired record gives 0, including records with neither class recorded. Records with missing counts in either class are reported separately.',
    'Recorded overlap is 100 times the joint-detection count divided by the paired denominator. Forest contrasts are Mbololo minus Ngangao in percentage points. A zero paired denominator produces an undefined percentage.',
    'The comparison is descriptive for the surveyed neighbourhoods. The available sampling documentation establishes neither record independence nor a cluster structure suitable for resampling, so confidence intervals are omitted. Wilson and Newcombe intervals require a defensible independent-record design.',
    'Recorded detections describe the survey evidence. The analysis does not estimate detection-corrected occupancy or test whether the size classes overlap more often than expected from their marginal frequencies.',
    '',
    'Independent-proportion intervals: Newcombe (1998), https://doi.org/10.1002/%28sici%291097-0258%2819980430%2917%3A8%3C873%3A%3Aaid-sim779%3E3.0.co%3B2-i',
    'Detection and occupancy: MacKenzie et al. (2002), https://armi.usgs.gov/content/?contentid=190',
    '',
    'Figure 9 shows the four-category breakdown. Table 4 reports the four size-class categories. Additional overlap estimates and forest contrasts are unnumbered diagnostics.')
  if (!provisional) writeLines(methods, output('METHODS.md'))
  results <- c('# Results. Recorded size-class overlap', '')
  for (i in seq_len(nrow(overlap))) {
    d <- overlap[i,]
    percent <- if (is.finite(d$recorded_overlap_percent))sprintf('%.1f%%', d$recorded_overlap_percent)else 'undefined'
    results <- c(results, sprintf('- %s, %s: %d/%d (%s); %d records missing either size class.',
        d$species, d$forest, d$both_size_classes_records, d$paired_eligible_records,
        percent, d$missing_either_class_records))
  }
  for (i in seq_len(nrow(contrast))) {
    d <- contrast[i,]
    difference <- d$difference_mbololo_minus_ngangao_pp
    results <- c(results, if (is.finite(difference))
      sprintf('- %s: Mbololo minus Ngangao = %+.1f percentage points.', d$species, difference)else
      sprintf('- %s: the difference is undefined because one or both fragments have no paired records.', d$species))
  }
  if (!provisional) writeLines(results, output('RESULTS.md'))
  invisible(list(overlap = overlap, contrast = contrast))
}

# ============================================================
# 5. Jao provisional
# ============================================================

jao_provisional <- function(root, here, out, absent) {
  path <- file.path(root, 'analysis/Datasets/derived/worksheet_focal_count_ledger.csv')
  if (!file.exists(path)) {
    env <- new.env(parent = environment())
    source(file.path(here, '03_focal_abundance.R'), local = env)
    env$run_population_distribution(root, here)
  }
  if (!file.exists(path))stop('Question 3 provisional count ledger unavailable: ', path, call. = FALSE)
  ledger <- read.csv(path, check.names = FALSE)
  key <- function(x) paste(x$record_id, x$species, sep = '|')
  j <- ledger[ledger$stage == 'Recorded juveniles',]
  l <- ledger[ledger$stage == 'Listed larger stems',]
  at <- match(key(j), key(l))
  stopifnot(nrow(j)>0, !anyDuplicated(key(j)), !anyDuplicated(key(l)), setequal(key(j), key(l)),
    identical(j$included, l$included[at]), identical(j$forest, l$forest[at]))
  z <- j[, c('record_id', 'forest', 'workbook', 'sheet', 'species', 'included', 'reason')]
  z$juvenile_count <- ifelse(j$stage_table_available, j$count, NA_real_)
  z$larger_count <- ifelse(l$stage_table_available[at], l$count[at], NA_real_)
  for (v in list(z$juvenile_count, z$larger_count))
  stopifnot(all(v[is.finite(v)] >= 0 & v[is.finite(v)] == floor(v[is.finite(v)])))
  z$cooccurrence <- jao_category(z$juvenile_count, z$larger_count)
  report <- jao_overlap_report(z[z$included,], out, provisional = TRUE)
  if ('possible_header_count' %in% names(l)) {
    alternative <- z[z$included,]
    alternative$larger_count <- l$possible_header_count[at][z$included]
    alternative$cooccurrence <- jao_category(alternative$juvenile_count, alternative$larger_count)
    groups <- split(alternative, list(alternative$forest, alternative$species), drop = TRUE)
    sensitivity <- do.call(rbind, lapply(groups, function(d) {
          paired <- sum(!is.na(d$cooccurrence));
          both <- sum(d$cooccurrence == 'Both stages', na.rm = TRUE)
          data.frame(forest = d$forest[1], species = d$species[1], paired_records = paired,
            both_records = both, recorded_overlap_percent = if (paired)100 * both / paired else NA_real_,
            scenario = 'Hypothetical unmatched focal-header inclusion')
        }))
    write.csv(sensitivity, rq_file(out, 'possible_header_overlap_sensitivity.csv'), row.names = FALSE, na = 'NA')
  }
  write.csv(z, rq_file(out, 'worksheet_record_stage_ledger.csv'), row.names = FALSE, na = 'NA')
  status <- c('Provisional Question 4 completed. Audited outputs remain unavailable.',
    paste('Count source:', path),
    'Recorded overlap uses all paired records, including neither-size-class records. Missing counts are reported separately.',
    'Counts use exact recorded labels; focal headers were not added.',
    'Forest differences are Mbololo minus Ngangao, in percentage points.',
    'Missing audited inputs:', absent)
  message(paste(status, collapse = '\n'))
  invisible(list(records = z, overlap = report$overlap, contrast = report$contrast, output = out, provisional = TRUE))
}

# ============================================================
# 6. Run juvenile abundance occurrence
# ============================================================

run_juvenile_abundance_occurrence <- function(root, refresh = FALSE) {
  root <- jao_project_root(root);
  here <- file.path(root, 'analysis/scripts')
  if (refresh) {
    env <- new.env(parent = environment())
    source(file.path(here, '03_focal_abundance.R'), local = env)
    env$run_population_distribution(root, here, refresh = FALSE)
  }
  src <- file.path(root, 'analysis/outputs/research question 3')
  out <- rq_dirs(root, 4)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  required <- rq_input(root, paste0('analysis/population_distribution/',
      c('focal_count_ledger', 'focal_header_evidence', 'taxonomic_decisions'), '.csv'))
  absent <- required[!file.exists(required)]
  if (length(absent))return(jao_provisional(root, here, out, absent))
  read <- function(n)read.csv(rq_input(root, paste0('analysis/population_distribution/', n, '.csv')), check.names = FALSE)
  put <- function(z, n)write.csv(z, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  ledger <- read('focal_count_ledger');
  headers <- read('focal_header_evidence')
  stopifnot(!anyDuplicated(ledger[c('record_id', 'species', 'stage')]))
  j <- ledger[ledger$stage == 'Recorded juveniles',]
  l <- ledger[ledger$stage == 'Listed larger stems',]
  key <- function(z)paste(z$record_id, z$species, sep = '|')
  at <- match(key(j), key(l));
  stopifnot(!anyNA(at))
  z <- j[, c('record_id', 'forest', 'source_file', 'worksheet', 'species', 'context', 'cultivated',
      'primary_record', 'exclusion_reason', 'equal_stage_search_area_verified', 'field_identity_verified')]
  z$juvenile_table_available <- j$stage_table_available
  z$larger_table_available <- l$stage_table_available[at]
  z$juvenile_count <- ifelse(z$juvenile_table_available, j$count, NA_real_)
  z$larger_count <- ifelse(z$larger_table_available, l$count[at], NA_real_)
  for (n in c('juvenile_count', 'larger_count')) {
    v <- z[[n]];
  stopifnot(all(v[is.finite(v)] >= 0), all(v[is.finite(v)] == floor(v[is.finite(v)])))
  }
  z$cooccurrence <- jao_category(z$juvenile_count, z$larger_count)
  z$paired_eligible <- !is.na(z$cooccurrence)
  z$missing_stage <- ifelse(z$paired_eligible, 'None', ifelse(!is.finite(z$juvenile_count) & !is.finite(z$larger_count),
      'Both', ifelse(!is.finite(z$juvenile_count), 'Juvenile', 'Larger stem')))
  put(z, 'record_stage_ledger');
  put(z[z$cultivated,], 'cultivated_records_excluded')
  put(z[!z$primary_record,], 'nonprimary_record_audit')
  put(z[!z$paired_eligible,], 'missing_stage_evidence')
  put(headers, 'focal_header_evidence_separate');
  put(read('taxonomic_decisions'), 'taxonomic_decisions')
  cohorts <- list(primary = z[z$primary_record & !z$cultivated,],
    broader_forest_context = z[!z$cultivated,])
  abund <- categories <- list();
  k <- 0L
  for (cohort in names(cohorts))for (forest in c('Ngangao', 'Mbololo'))for (sp in unique(z$species)) {
    v <- cohorts[[cohort]];
  v <- v[v$forest == forest & v$species == sp,]
    s <- jao_summarise(v);
    k <- k + 1L
    ids <- data.frame(cohort, forest, species = sp)
    abund[[k]] <- cbind(ids, s$abundance)
    categories[[k]] <- cbind(ids[rep(1, 4),], s$categories)
  }
  abund <- do.call(rbind, abund);
  categories <- do.call(rbind, categories)
  put(abund, 'juvenile_abundance_and_frequency');
  put(categories, 'stage_cooccurrence')
  primary <- abund[abund$cohort == 'primary',];
  put(primary, 'primary_summary')
  # Explicit changes from primary to broader; no ecological uncertainty interval implied.
  metrics <- c('source_records', 'eligible_juvenile_records', 'juvenile_total', 'juvenile_positive_records',
    'record_frequency_percent', 'paired_eligible_records', 'both_stages_records', 'larger_presence_given_juvenile')
  p <- abund[abund$cohort == 'primary', c('forest', 'species', metrics)]
  b <- abund[abund$cohort == 'broader_forest_context', c('forest', 'species', metrics)]
  delta <- merge(p, b, by = c('forest', 'species'), suffixes = c('_primary', '_broader'))
  for (m in metrics)delta[[paste0(m, '_change')]] <- delta[[paste0(m, '_broader')]] - delta[[paste0(m, '_primary')]]
  put(delta, 'broader_minus_primary_changes')
  # Hypothetical header additions affect Coffea larger-stem presence only.
  scenario <- z[!z$cultivated,]
  h <- headers$possible_unlisted_coffea_header[match(scenario$record_id, headers$plot_id)]
  stopifnot(!anyNA(h))
  scenario$hypothetical_header_addition <- as.integer(scenario$species == 'coffea fadenii' & h)
  scenario$larger_count <- scenario$larger_count + scenario$hypothetical_header_addition
  scenario$cooccurrence <- jao_category(scenario$juvenile_count, scenario$larger_count)
  scenario$scenario <- 'Unverified possible-header addition; not observed listed counts'
  put(scenario, 'possible_header_scenario_ledger')
  hs <- list();
  k <- 0L
  for (cohort in names(cohorts))for (forest in c('Ngangao', 'Mbololo'))for (sp in unique(z$species)) {
    v <- scenario[scenario$forest == forest & scenario$species == sp &
      (cohort != 'primary' | scenario$primary_record),]
    s <- jao_summarise(v);
    k <- k + 1L
    hs[[k]] <- cbind(data.frame(cohort, forest, species = sp, scenario = 'Unverified possible-header addition'), s$categories)
  }
  put(do.call(rbind, hs), 'possible_header_scenario_cooccurrence')
  report <- jao_overlap_report(cohorts$primary, out)
  jao_reporting(out)
  inputs <- c(file.path(here, '04_size_class_overlap.R'),
    rq_input(root, paste0('analysis/population_distribution/', c('focal_count_ledger.csv',
          'focal_header_evidence.csv', 'taxonomic_decisions.csv', 'input_and_code_hashes.csv'))))
  put(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))), 'input_and_code_hashes')
  writeLines(capture.output(sessionInfo()), rq_file(out, 'sessionInfo.txt'))
  message('Juvenile abundance and occurrence complete: ', out)
  invisible(list(records = z, abundance = abund, categories = categories, overlap = report$overlap, contrast = report$contrast, output = out))
}

# ============================================================
# 7. Jao reporting
# ============================================================

# Publication naming continues Figures 1-3, Tables 1-3, Figures S1-S7 and Tables S1-S10.
jao_reporting <- function(out) {
  read <- function(n)read.csv(rq_file(out, paste0(n, '.csv')), check.names = FALSE)
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  a <- read('juvenile_abundance_and_frequency');
  c <- read('stage_cooccurrence')
  cats <- c('Both stages', 'Juveniles only', 'Larger stems only', 'Neither stage listed')
  keys <- c('both_stages', 'juveniles_only', 'larger_stems_only', 'neither_stage_listed')
  for (i in seq_along(cats)) {
    cc <- c[c$category == cats[i],]
    at <- match(paste(a$cohort, a$forest, a$species), paste(cc$cohort, cc$forest, cc$species))
    stopifnot(!anyNA(at))
    a[[paste0(keys[i], '_n')]] <- cc$records[at]
    a[[paste0(keys[i], '_percent')]] <- cc$percent[at]
  }
  put(read('record_stage_ledger'), 'Table_S11_record_level_juvenile_and_larger_stem_evidence')
  put(a, 'Table_S12_primary_and_broader_cohort_stage_summaries')
  h <- read('possible_header_scenario_cooccurrence')
  at <- match(paste(h$cohort, h$forest, h$species, h$category), paste(c$cohort, c$forest, c$species, c$category))
  stopifnot(!anyNA(at))
  h$listed_records <- c$records[at];
  h$listed_percent <- c$percent[at]
  h$change_in_records <- h$records - h$listed_records;
  h$change_in_percentage_points <- h$percent - h$listed_percent
  put(h, 'Table_S13_possible_header_cooccurrence_sensitivity')
}

# ============================================================
# 8. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep('^--file=', commandArgs(), value = TRUE)
  start <- if (length(arg))dirname(normalizePath(sub('^--file=', '', arg[1]), mustWork = TRUE))else getwd()
  root <- jao_project_root(start)
  run_juvenile_abundance_occurrence(root)
  source(file.path(file.path(root, 'analysis'), 'scripts', 'combine_output_tables.R'), local = TRUE)
  combine_output_tables(file.path(file.path(root, 'analysis'), 'outputs', 'research question 4'), '04_size_class_overlap.R')
}
