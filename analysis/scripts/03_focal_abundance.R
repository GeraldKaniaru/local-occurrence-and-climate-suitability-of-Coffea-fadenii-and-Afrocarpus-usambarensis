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
# 1. Pd project root
# ============================================================

# Research question 3: For each focal species, how do recorded abundance, occurrence frequency, and concentration among survey records vary by fragment and size class?
# Research question 3: Recorded abundance, frequency and concentration among Taita source records.
pd_project_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = '/', mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, 'analysis/scripts/03_focal_abundance.R')))return(p)
    parent <- dirname(p);
    if (identical(parent, p))break;
    p <- parent
  }
  stop('Set the working directory to the Taita project or a folder within it.', call. = FALSE)
}

# ============================================================
# 2. Pd temp ready
# ============================================================

pd_temp_ready <- function() {
  tmp <- tempdir(check = TRUE)
  if (!dir.exists(tmp))dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  probe <- tempfile('population_write_check_', tmpdir = tmp)
  ok <- suppressWarnings(file.create(probe));
  if (isTRUE(ok))unlink(probe)
  if (!isTRUE(ok))stop('R temporary directory is unavailable. Restart R and Source the script again.', call. = FALSE)
  invisible(tmp)
}

# ============================================================
# 3. Population reporting
# ============================================================

pd_table_3 <- function(tab) {
  tab <- tab[order(match(tab$species, c('coffea fadenii', 'afrocarpus usambarensis')),
    match(tab$stage, c('Recorded juveniles', 'Listed larger stems')),
    match(tab$forest, c('Ngangao', 'Mbololo'))), ]
  format_value <- function(x, digits) ifelse(is.finite(x), sprintf('%.*f', digits, x), '\u2014')
  data.frame(
    'Species' = paste0(toupper(substr(tab$species, 1, 1)), '. ',
      sub('^[^ ]+ ', '', tab$species)),
    'Size class' = ifelse(tab$stage == 'Recorded juveniles', 'Juveniles', 'Listed larger stems'),
    'Forest' = tab$forest,
    'Recorded count' = tab$total_stems,
    'Frequency (%)' = format_value(tab$record_frequency_percent, 1),
    'Mean count across all records' = format_value(tab$mean_per_eligible_record, 2),
    'Mean count among records with detections' = format_value(tab$mean_per_positive_record, 2),
    check.names = FALSE)
}

pd_write_table_3 <- function(table, path) {
  lines <- capture.output(write.csv(table, row.names = FALSE, na = 'NA'))
  lines <- gsub('"<U+2014>"', '"\u2014"', lines, fixed = TRUE)
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}

population_reporting <- function(out, primary, summary, ledger, headers, location) {
  table_3 <- pd_table_3(transform(primary, total_stems = total_listed_count,
    mean_per_eligible_record = mean_per_record,
    mean_per_positive_record = ifelse(records_with_listed_count > 0,
      total_listed_count / records_with_listed_count, NA_real_)))
  pd_write_table_3(table_3, rq_file(out, 'Table_1_Abundance_and_occurrence_by_fragment_and_size_class.csv'))
  tables <- list(Table_S7_focal_counts_and_record_eligibility = ledger,
    Table_S8_cohort_sensitivity_and_concentration = summary,
    Table_S9_focal_header_evidence = headers,
    Table_S10_worksheet_location_audit = location)
  for (n in names(tables))write.csv(tables[[n]], rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  fmt <- function(x)ifelse(is.finite(x), formatC(x, digits = 2, format = 'f'), 'NA')
  lines <- c('# Table 3. Recorded focal abundance and frequency', '',
    '| Species | Forest | Stage | Eligible records | Total | Median [Q1, Q3] | Range | Positive/eligible (%) |',
    '|---|---|---|---:|---:|---|---|---|')
  for (i in seq_len(nrow(primary))) {z <- primary[i,];
    vals <- c(z$species, z$forest, z$stage, z$eligible_records, z$total_listed_count,
      paste0(fmt(z$median_per_record), ' [', fmt(z$q1), ', ', fmt(z$q3), ']'), paste0(z$minimum, '-', z$maximum),
      paste0(z$records_with_listed_count, '/', z$eligible_records, ' (', fmt(z$record_frequency_percent), '%)'))
    lines <- c(lines, paste0('| ', paste(vals, collapse = ' | '), ' |'))}
  writeLines(c(lines, '', 'Counts describe selected survey records, not population density or detection-corrected occupancy.'), rq_file(out, 'Table_1_Abundance_and_occurrence_by_fragment_and_size_class.md'))
  figures <- c('Figure_3_recorded_focal_abundance', 'Figure_S5_recorded_occurrence_frequency',
    'Figure_S6_concentration_among_survey_records', 'Figure_S7_worksheet_source_record_locations')
  register <- data.frame(label = c('Figure 3', 'Figure S5', 'Figure S6', 'Figure S7', 'Table 3', 'Table S7', 'Table S8', 'Table S9', 'Table S10'), file = rq_relative(c(paste0(figures, '.pdf'), paste0(names(tables), '.csv'))))
  write.csv(register, rq_file(out, 'publication_output_register.csv'), row.names = FALSE)
  writeLines(c('# Population-distribution reporting', '',
      'Figure 3: Individual-record listed focal counts, separately by species and operational stage; forest assignment is indicated by colour.',
      'Figure S5: Record frequencies, with explicit positive-record numerators and eligible-record denominators.',
      'Figure S6: Cumulative shares of listed counts among records ranked separately within each stage. These are not geographic clustering or habitat-loss redundancy curves.',
      'Figure S7: Worksheet coordinates; forest, primary eligibility and cultivated contexts are distinguished. Positional accuracy and datum remain unconfirmed.',
      'Table 3: Primary abundance and record frequency.',
      'Table S7: Complete focal count ledger and record eligibility.',
      'Table S8: Cohort sensitivity, totals, observed and count-standardised concentration.',
      'Table S9: Separate focal-header evidence; hypothetical additions do not enter listed counts.',
      'Table S10: Coordinate provenance and location audit.'), rq_file(out, 'REPORTING_GUIDE.md'))
}

# ============================================================
# 4. Pd concentration
# ============================================================

pd_concentration <- function(x) {
  x <- x[is.finite(x)]
  if (any(x<0 | x != floor(x)))stop('Counts must be nonnegative integers')
  n <- if (length(x))sum(x)else NA_real_
  positive <- is.finite(n) && n>0
  c(total = n, eligible_records = length(x), records_with_count = sum(x>0),
    concentration = if (positive)sum((x / n)^2)else NA_real_,
    effective_number_of_survey_records = if (positive)1 / sum((x / n)^2)else NA_real_,
    largest_record_share = if (positive)max(x) / n else NA_real_,
    distinct_pair_concentration = if (positive && n >= 2)sum(x * (x - 1)) / (n * (n - 1))else NA_real_)
}

# ============================================================
# 5. Pd equal count expectation
# ============================================================

pd_equal_count_expectation <- function(x, m) {
  x <- x[is.finite(x)];
  n <- sum(x)
  if (!is.finite(m) || m<1 || m != floor(m) || m>n)return(NA_real_)
  if (m == 1)return(1)
  # Exact expected empirical concentration after a simple random subsample of
  # m observed individuals WITHOUT replacement. No ecological null is assumed.
  c <- pd_concentration(x)['distinct_pair_concentration']
  unname(1 / m + (m - 1) * c / m)
}

# ============================================================
# 6. Pd two part report
# ============================================================

pd_two_part_report <- function(ledger, selected, out, provisional = FALSE) {
  forests <- c('Ngangao', 'Mbololo')
  species <- c('coffea fadenii', 'afrocarpus usambarensis')
  stages <- c('Recorded juveniles', 'Listed larger stems')
  stage_labels <- c('Juveniles', 'Larger stems')
  colours <- c('#D4B58A', '#508A92')
  symbols <- c(16, 17)
  z <- ledger[which(selected & ledger$forest %in% forests),]
  if (anyDuplicated(z[c('record_id', 'species', 'stage')]))
  stop('Duplicate record-species-stage count in Question 3')
  eligible <- z$stage_table_available %in% TRUE & is.finite(z$count)
  if (any(z$count[eligible]<0 | z$count[eligible] != floor(z$count[eligible])))
  stop('Question 3 requires nonnegative integer counts')
  z$count[!eligible] <- NA_real_
  prefix <- if (provisional)'provisional_'else ''
  output <- function(name)rq_file(out, if (grepl('[.]md$', name))paste0(prefix, name)else name)
  put <- function(x, name)write.csv(x, output(paste0(name, '.csv')), row.names = FALSE, na = 'NA')
  groups <- expand.grid(forest = forests, stage = stages, species = species, stringsAsFactors = FALSE)
  summary <- vector('list', nrow(groups))
  ranked <- curves <- vector('list', nrow(groups))
  for (i in seq_len(nrow(groups))) {
    g <- groups[i,]
    d <- z[z$forest == g$forest & z$species == g$species & z$stage == g$stage,]
    d <- d[order(d$count, d$record_id, na.last = NA),]
    p <- d[d$count>0,]
    n <- nrow(d)
    k <- nrow(p)
    total <- sum(p$count)
    gini <- if (k >= 3)sum((2 * seq_len(k) - k - 1) * p$count) / (k * total)else NA_real_
    summary[[i]] <- data.frame(g, selected_records = sum(z$forest == g$forest &
        z$species == g$species & z$stage == g$stage), eligible_records = n,
      zero_records = n - k, positive_records = k, total_stems = if (n)total else NA_real_,
      mean_per_eligible_record = if (n)total / n else NA_real_,
      record_frequency = if (n)k / n else NA_real_,
      record_frequency_percent = if (n)100 * k / n else NA_real_,
      mean_per_positive_record = if (k)total / k else NA_real_, gini_positive = gini,
      concentration_status = if (k >= 3)'Lorenz and empirical Gini'else
      if (k)'Sparse; show individual counts'else 'No positive counts')
    ranked[[i]] <- data.frame(p[c('forest', 'stage', 'species', 'record_id', 'count')],
      rank = seq_len(k))
    if (k >= 3)curves[[i]] <- data.frame(g[rep(1, k + 1),], positive_records = k,
      cumulative_record_share = (0:k) / k, cumulative_stem_share = c(0, cumsum(p$count) / total))
  }
  tab <- do.call(rbind, summary)
  tab$unavailable_records <- tab$selected_records - tab$eligible_records
  contrasts <- merge(tab[tab$forest == forests[1],], tab[tab$forest == forests[2],],
    by = c('species', 'stage'), suffixes = c('_ngangao', '_mbololo'))
  delta <- contrasts[c('species', 'stage')]
  for (metric in c('total_stems', 'mean_per_eligible_record', 'record_frequency_percent',
      'mean_per_positive_record', 'gini_positive'))
  delta[[paste0(metric, '_mbololo_minus_ngangao')]] <-
  contrasts[[paste0(metric, '_mbololo')]] - contrasts[[paste0(metric, '_ngangao')]]
  rank_table <- do.call(rbind, ranked)
  curve_table <- if (any(lengths(curves)>0L))do.call(rbind, curves)else
  data.frame(groups[FALSE,], positive_records = integer(),
    cumulative_record_share = numeric(), cumulative_stem_share = numeric())
  pd_write_table_3(pd_table_3(tab), output('Table_1_Abundance_and_occurrence_by_fragment_and_size_class.csv'))
  put(delta, 'forest_contrasts')
  put(rank_table, 'positive_record_ranked_counts')
  put(curve_table, 'positive_record_lorenz_coordinates')
  put(z[is.na(z$count),], 'unavailable_focal_counts')
  species_label <- function(sp) {
    name <- paste0(toupper(substr(sp, 1, 1)), substring(sp, 2))
    bquote(italic(.(name)))
  }
  draw_abundance <- function() {
    panels <- list()
    for (sp in species) {
      d <- z[z$species == sp & is.finite(z$count), ]
      t <- tab[tab$species == sp, ]
      t$x <- match(t$stage, stages) + c(.17, -.17)[match(t$forest, forests)]
      d$x <- match(d$stage, stages) + c(.17, -.17)[match(d$forest, forests)]
      if (nrow(d)) {
        d$x <- d$x + ave(seq_len(nrow(d)), d$forest, d$stage, d$count,
          FUN = function(ids) if (length(ids) == 1L) 0 else seq(-.09, .09, length.out = length(ids)))
      }
      for (occurrence in c(FALSE, TRUE)) {
        index <- match(sp, species) + if (occurrence) 2L else 0L
        if (index == 4L) t$x <- 2 * match(t$stage, stages) - t$x
        subtitle <- bquote(.(paste0('(', letters[index], ')')) ~ .(species_label(sp)))
        minimum <- 0
        maximum <- if (index == 1L) max(30, 20 * ceiling(max(c(0, d$count)) / 20)) else
          if (occurrence) 100 else max(c(1, d$count)) * 1.12
        breaks <- if (index == 1L) seq(0, maximum, if (maximum > 30) 20 else 5) else
          if (occurrence) seq(0, 100, 20) else ggplot2::waiver()
        plot <- ggplot2::ggplot() +
          ggplot2::scale_x_continuous(limits = c(.55, 2.45), breaks = 1:2, labels = stage_labels) +
          ggplot2::scale_y_continuous(breaks = breaks,
            expand = ggplot2::expansion(mult = c(.035, .075))) +
          ggplot2::coord_cartesian(ylim = c(minimum, maximum)) +
          ggplot2::annotate('text', x = .57, y = maximum + .06 * (maximum - minimum), label = as.expression(subtitle),
            hjust = 0, vjust = 1, size = 6.2 - 4 * 25.4/72) +
          ggplot2::scale_fill_manual(values = setNames(colours, forests), breaks = forests) +
          ggplot2::scale_shape_manual(values = setNames(c(21, 24), forests), breaks = forests) +
          ggplot2::labs(x = NULL, y = if (occurrence) 'Records with detections (%)' else
            'Stems per eligible record', fill = NULL, shape = NULL) +
          ggplot2::theme_bw(base_size = 13) +
          ggplot2::theme(text = ggplot2::element_text(colour = 'black'),
            axis.text = ggplot2::element_text(size = 13, colour = 'black'),
            axis.text.x = ggplot2::element_text(face = 'bold'),
            axis.title.y = ggplot2::element_text(size = 14, face = 'bold', margin = ggplot2::margin(r = 10)),
            legend.position = 'bottom', legend.direction = 'horizontal',
            legend.text = ggplot2::element_text(size = 14),
            legend.key.size = grid::unit(8, 'mm'),
            legend.background = ggplot2::element_rect(fill = 'white', colour = NA))
        if (occurrence) {
          t$label <- ifelse(t$eligible_records > 0, sprintf('%.1f%%', t$record_frequency_percent), '0 eligible')
          label_gap <- .055 * maximum
          t$label_y <- ifelse(t$eligible_records > 0, t$record_frequency_percent + label_gap, label_gap)
          if (index == 3L) {
            high <- which(t$record_frequency_percent >= 90)
            t$label_y[high] <- t$record_frequency_percent[high] - label_gap
          }
          plot <- plot + ggplot2::geom_rect(data = t[which(t$record_frequency_percent >= minimum), ],
            ggplot2::aes(xmin = x - .125, xmax = x + .125, ymin = minimum,
              ymax = record_frequency_percent, fill = forest),
            colour = 'black', linewidth = .45) +
            ggplot2::geom_text(data = t, ggplot2::aes(x, label_y, label = label), size = 5.3 - 2 * 25.4/72)
        } else {
          plot <- plot + ggplot2::geom_point(data = d,
            ggplot2::aes(x, count, fill = forest, shape = forest), size = 3.6) +
            ggplot2::geom_segment(data = t[is.finite(t$mean_per_eligible_record), ],
              ggplot2::aes(x = x - .149, xend = x + .149, y = mean_per_eligible_record,
                yend = mean_per_eligible_record, colour = forest), linewidth = 2.33,
              lineend = 'square') +
            ggplot2::scale_colour_manual(values = setNames(colours, forests), guide = 'none')
        }
        if (index <= 2L) {
          plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(5.5, 5.5, 1 + 2 * 72/25.4, 5.5))
        } else {
          plot <- plot + ggplot2::theme(plot.margin = ggplot2::margin(1, 5.5, 1, 5.5))
        }
        if (index == 1L) plot <- plot + ggplot2::guides(
          shape = ggplot2::guide_legend(override.aes = list(size = 5.2)))
        if (index %in% c(2L, 4L)) plot <- plot + ggplot2::labs(y = NULL)
        if (index != 1L) plot <- plot + ggplot2::guides(fill = 'none', shape = 'none')
        panels[[index]] <- plot
      }
    }
    print(patchwork::wrap_plots(panels, ncol = 2, guides = 'collect') &
      ggplot2::theme(legend.position = 'bottom',
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.margin = ggplot2::margin(2, 0, 0, 0)))
  }
  draw_concentration <- function() {
    forest_order <- c('Mbololo', 'Ngangao')
    palette <- setNames(colours, forests)
    focal <- species[1]
    counts <- tab[tab$species == focal, ]
    panels <- lapply(seq_along(stages), function(j) {
      d <- counts[counts$stage == stages[j] & is.finite(counts$gini_positive), ]
      d$label_y <- c(.75, .63)[match(d$forest, forest_order)]
      curves <- curve_table[curve_table$species == focal & curve_table$stage == stages[j], ]
      heading <- bquote(.(paste0('(', letters[j], ')')) ~ .(species_label(focal)) ~ '|' ~ bold(.(stage_labels[j])))
      ggplot2::ggplot() +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey55') +
        ggplot2::geom_line(data = curves, ggplot2::aes(cumulative_record_share,
          cumulative_stem_share, colour = forest), linewidth = 1.7) +
        ggplot2::annotate('text', x = .025, y = .975, label = as.expression(heading),
          hjust = 0, vjust = 1, size = 5.4) +
        ggplot2::geom_label(data = d, ggplot2::aes(x = .06, y = label_y, fill = forest,
          label = sprintf('Gini = %.2f', gini_positive)), hjust = 0, size = 5.3,
          colour = 'black', fontface = 'bold', label.r = grid::unit(2, 'mm'), show.legend = FALSE) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = 'Cumulative share of positive records',
          y = 'Cumulative share of stems', colour = NULL) +
        ggplot2::scale_colour_manual(values = palette, limits = forest_order) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linewidth = 5.1))) +
        ggplot2::scale_fill_manual(values = palette) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(text = ggplot2::element_text(colour = 'black'),
          axis.text = ggplot2::element_text(size = 16, colour = 'black'),
          axis.title = ggplot2::element_text(size = 18, face = 'bold'),
          legend.position = 'bottom', legend.direction = 'horizontal',
          legend.text = ggplot2::element_text(size = 17),
          plot.margin = ggplot2::margin(5.5, 5.5, 1, 5.5))
    })
    print(patchwork::wrap_plots(panels, ncol = 2, guides = 'collect') &
      ggplot2::theme(legend.position = 'bottom',
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.margin = ggplot2::margin(2, 0, 0, 0)))
  }
  plots <- list(Figure_5_Abundance_and_occurrence_of_focal_species = draw_abundance,
    Figure_S8_Lorenz_curves_for_Coffea_fadenii = draw_concentration)
  for (name in names(plots))for (ext in c('pdf', 'png')) {
    path <- output(paste0(name, '.', ext))
    width <- if (name == 'Figure_5_Abundance_and_occurrence_of_focal_species') 10.5 else 13.5
    height <- if (name == 'Figure_5_Abundance_and_occurrence_of_focal_species') 6.5 else 6
    if (ext == 'pdf') pdf(path, width = width, height = height, compress = dir.exists(tempdir()))
    else png(path, width = width * 220, height = height * 220, res = 220)
    tryCatch(plots[[name]](), finally = dev.off())
  }
  if (!provisional) writeLines(c('# Question 3 methods and figure guide', '',
      if (provisional)'Provisional worksheet counts from the primary forest-associated selection. Only exact focal labels after case and whitespace normalisation enter counts. Focal headers and cultivated records are outside the primary comparison.'else
      'Primary audited forest-associated records. Focal headers and cultivated records are outside the primary comparison.',
      '',
      'Counts are defined for each record, species and operational size class. A usable stage table with no listed focal count contributes zero. Unavailable or unusable count evidence contributes NA and stays outside the denominator.',
      'Table 3 reports totals, eligible and positive record counts, occurrence frequency and means with and without eligible zeros. For groups with detections, mean per eligible record equals record frequency multiplied by mean per positive record. Conditional means are undefined when there are no detections.',
      'Figure 7 shows individual counts, horizontal mean markers and occurrence percentages. Horizontal offsets separate identical counts visually. Figure 8 shows Lorenz curves and empirical Gini coefficients for groups with at least three positive records. Sparse groups are shown as ranked dots when both fragments are sparse, or as individual counts in the legend when the other fragment has a curve. The positive-record ranking table gives record IDs.',
      'The Gini is the mean absolute pairwise count difference divided by twice the mean, using positive records only and no small-sample correction. It measures inequality in recorded counts among detections. Forest contrasts are Mbololo minus Ngangao; occurrence differences are in percentage points.',
      'Surveyed areas and the independence of sampling units are unresolved. Results are descriptive for the surveyed neighbourhoods, use counts per record and include no density extrapolation, detection correction, geographic clustering test or formal two-part model.',
      '',
      'Abundance and occurrence context: Gaston et al. (2000), https://doi.org/10.1046/j.1365-2664.2000.00485.x',
      'Lorenz and Gini precedent in plant data: Weiner and Solbrig (1984), https://doi.org/10.1007/BF00379630'),
    output('QUESTION_3_METHODS.md'))
  invisible(tab)
}

# ============================================================
# 7. Pd provisional
# ============================================================

pd_provisional <- function(root, here, out, absent) {
  env <- new.env(parent = environment())
  source(file.path(here, '01_woody_composition.R'), local = env)
  input <- env$floristic_provisional(root, out, export = FALSE)
  rec <- input$records
  rows <- input$source_rows
  species <- c('coffea fadenii', 'afrocarpus usambarensis')
  stages <- c('Listed larger stems', 'Recorded juveniles')
  ledger <- do.call(rbind, lapply(seq_len(nrow(rec)), function(i) {
        r <- rec[i,]
        do.call(rbind, lapply(stages, function(stage) {
              z <- rows[rows$record_id == r$record_id & rows$stage == stage,]
              ev <- input$evidence[input$evidence$record_id == r$record_id & input$evidence$stage == stage,]
              stopifnot(nrow(ev) == 1L)
              do.call(rbind, lapply(species, function(sp) {
                    x <- z[z$recorded_label == sp,]
                    available <- ev$status %in% c('Recorded counts available', 'Surveyed: explicit zero counts')
                    count <- if (available)sum(x$count)else NA_real_
                    numeric_label <- function(v) {
                      v <- sub('[[:space:]]*cm$', '', trimws(v), ignore.case = TRUE)
                      suppressWarnings(as.numeric(sub(',', '.', v, fixed = TRUE)))
                    }
                    focal <- isTRUE(r$focal_header_taxon == sp)
                    header <- if (focal)numeric_label(if (stage == stages[1])r$focal_header_dbh else r$focal_header_saplings)else NA_real_
                    dbh_matches <- if (stage == stages[1] && is.finite(header) && header>0)
                    sum(abs(numeric_label(x$raw_dbh) - header)<1e-8, na.rm = TRUE)else NA_integer_
                    possible_addition <- focal && stage == stages[1] && available && is.finite(header) && header>0 && dbh_matches == 0L
                    data.frame(r, species = sp, stage, count, stage_table_available = available,
                      stage_table_status = ev$status, matching_rows = nrow(x), invalid_matching_rows = sum(!x$usable_count),
                      focal_header_value = header, listed_same_taxon_dbh_matches = dbh_matches,
                      possible_header_addition = possible_addition,
                      possible_header_count = if (available)count + as.integer(possible_addition)else NA_real_,
                      juvenile_header_minus_list = if (stage == stages[2] && is.finite(header) && available)header - count else NA_real_)
                  }))
            }))
      }))
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  put(rec, 'worksheet_selection');
  put(ledger, 'worksheet_focal_count_ledger')
  derived <- file.path(root, 'analysis/Datasets/derived')
  dir.create(derived, recursive = TRUE, showWarnings = FALSE)
  write.csv(ledger, file.path(derived, 'worksheet_focal_count_ledger.csv'), row.names = FALSE, na = 'NA')
  summary <- pd_two_part_report(ledger, ledger$included, out, provisional = TRUE)
  focal <- ledger[which(ledger$species == ledger$focal_header_taxon),]
  put(focal, 'focal_header_vs_list_audit')
  sensitivity <- do.call(rbind, lapply(split(ledger[ledger$included,],
        list(ledger$forest[ledger$included], ledger$species[ledger$included], ledger$stage[ledger$included]), drop = TRUE), function(z)
      data.frame(forest = z$forest[1], species = z$species[1], stage = z$stage[1],
        eligible_records = sum(z$stage_table_available), listed_total = sum(z$count, na.rm = TRUE),
        unmatched_header_candidates = sum(z$possible_header_addition),
        possible_header_total = sum(z$possible_header_count, na.rm = TRUE))))
  put(sensitivity, 'focal_header_count_sensitivity')
  status <- c('Provisional Question 3 completed; audited inputs remain unavailable.',
    'Question 1 worksheet selection was applied. Only exact focal labels after case and whitespace normalisation were counted.',
    'Counts use numbered larger-stem and woody juvenile table rows. Focal headers were not added.',
    'A zero requires a usable stage table. Missing table evidence or any invalid numbered count yields NA.',
    'The header-versus-list audit reports juvenile count discrepancies and candidate same-name/same-DBH matches. The possible-header alternative adds at most one measured focal header per eligible record with no exact listed match. Its totals are hypothetical.',
    'Frequency uses records with available counts. Concentration uses records with detections; Gini is shown for three or more detections.',
    'Missing audited inputs:', absent)
  message(paste(status, collapse = '\n'))
  invisible(list(ledger = ledger, summary = summary, output = out, provisional = TRUE))
}

# ============================================================
# 8. Run population distribution
# ============================================================

run_population_distribution <- function(root, here = file.path(root, 'analysis/scripts'), refresh = FALSE) {
  root <- pd_project_root(root);
  here <- file.path(root, 'analysis/scripts');
  pd_temp_ready()
  local_lib <- file.path(here, '.r-library')
  if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
  needed <- c('readxl', 'dplyr', 'jsonlite')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))stop('Missing packages: ', paste(missing, collapse = ', '), '. Install explicitly.')
  # Reuse existing relevant analysis files; no config or support folder.
  source_paths <- rq_input(root, c('analysis/vegetation_structure/record_eligibility.csv',
      'analysis/vegetation_structure/individual_stem_measurement_audit.csv',
      'analysis/floristic_composition/juvenile_source_row_audit.csv',
      'analysis/vegetation_structure/taxonomic_decisions.csv'))
  absent <- source_paths[!file.exists(source_paths)]
  out <- rq_dirs(root, 3)
  if (length(absent))return(pd_provisional(root, here, out, absent))
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE, na = 'NA')
  q <- read.csv(source_paths[1], check.names = FALSE)
  a <- read.csv(source_paths[2], check.names = FALSE)
  j <- read.csv(source_paths[3], check.names = FALSE)
  tax <- read.csv(source_paths[4], check.names = FALSE)
  stopifnot(!anyDuplicated(q$plot_id), !anyNA(match(a$plot_id, q$plot_id)),
    !anyNA(match(j$plot_id, q$plot_id)))
  species <- c('coffea fadenii', 'afrocarpus usambarensis')
  stages <- c('Listed larger stems', 'Recorded juveniles')
  display <- function(x)paste0(toupper(substr(x, 1, 1)), substring(x, 2))
  ledger <- list();
  k <- 0L
  for (i in seq_len(nrow(q)))for (sp in species)for (stage in stages) {
    z <- if (stage == stages[1])a else j
    source <- z[z$plot_id == q$plot_id[i],, drop = FALSE]
    focal <- source[source$analysis_name == sp & source$woody_eligible,, drop = FALSE]
    available <- nrow(source)>0
    count <- if (!available)NA_real_ else if (stage == stages[1])nrow(focal)else sum(focal$count)
    if (is.finite(count) && (count<0 || count != floor(count)))stop('Invalid count for ', q$plot_id[i])
    k <- k + 1L;
    ledger[[k]] <- data.frame(record_id = q$plot_id[i], forest = q$forest[i],
      source_file = q$source_file[i], worksheet = q$sheet[i], species = sp, stage,
      context = q$cohort[i], cultivated = q$cultivated[i], primary_record = q$primary_eligible[i],
      exclusion_reason = q$exclusion_reason[i], stage_source_rows = nrow(source),
      stage_table_available = available, count = count, listed_occurrence = if (is.finite(count))count>0 else NA,
      count_status = if (!available)'Missing stage-table evidence; count unknown'else
      if (count == 0)'No matching listed focal count; not verified absence'else
      'Listed table count; individuals and detection not certified',
      coffea_header_count_disputed = if (sp == species[1])q$count_disputed[i]else NA,
      equal_stage_search_area_verified = q$equal_stage_search_area_verified[i],
      field_identity_verified = q$field_identity_verified[i])
  }
  ledger <- dplyr::bind_rows(ledger)
  put(ledger, 'focal_count_ledger');
  put(q, 'record_reconciliation')
  write.csv(ledger, rq_input(root, 'focal_count_ledger.csv'), row.names = FALSE, na = 'NA')
  put(tax, 'taxonomic_decisions')
  put(ledger[!ledger$primary_record,], 'excluded_and_cultivated_count_records')
  put(ledger[!ledger$stage_table_available,], 'missing_stage_table_evidence')
  # Header information remains a different evidence stream, never a count added
  # to listed plants. Same species and diameter indicates only a candidate match.
  headers <- q[, c('plot_id', 'forest', 'cohort', 'primary_eligible', 'cultivated',
      'focal_species_header', 'focal_dbh_header', 'juvenile_header', 'count_disputed',
      'matching_header_rows', 'header_duplicate_status', 'possible_unlisted_coffea_header',
      'coffea_possible_header_scenario')]
  headers$interpretation <- 'Separate header evidence; hypothetical additions never enter listed-count summaries'
  put(headers, 'focal_header_evidence')
  write.csv(headers, rq_input(root, 'focal_header_evidence.csv'), row.names = FALSE, na = 'NA')
  header_scenarios <- list();
  hi <- 0L
  for (cohort in c('primary', 'broader_forest_context'))for (f in c('Ngangao', 'Mbololo')) {
    zz <- q[q$forest == f & if (cohort == 'primary')q$primary_eligible else !q$cultivated,]
    hi <- hi + 1L;
    header_scenarios[[hi]] <- data.frame(cohort, forest = f, survey_records = nrow(zz),
      listed_coffea_stems = sum(zz$adult_cf), possible_header_coffea_stems = sum(zz$coffea_possible_header_scenario),
      hypothetical_additions = sum(zz$coffea_possible_header_scenario) - sum(zz$adult_cf),
      status = 'Unverified inclusion alternative; never substituted for observed table counts')
  }
  put(dplyr::bind_rows(header_scenarios), 'possible_header_sensitivity')
  cohorts <- list(primary = q$plot_id[q$primary_eligible],
    broader_forest_context = q$plot_id[!q$cultivated], cultivated = q$plot_id[q$cultivated])
  summary <- contributions <- list();
  si <- ci <- 0L
  for (cohort in names(cohorts))for (forest in c('Ngangao', 'Mbololo'))for (sp in species) {
    ids <- intersect(cohorts[[cohort]], q$plot_id[q$forest == forest])
    if (!length(ids))next
    z <- ledger[ledger$record_id %in% ids & ledger$species == sp,]
    # A fairer arithmetic comparison uses exactly the same records in both stages.
    az <- z[z$stage == stages[1],];
    bz <- z[z$stage == stages[2],]
    common <- intersect(az$record_id[is.finite(az$count)], bz$record_id[is.finite(bz$count)])
    av <- az$count[match(common, az$record_id)];
    bv <- bz$count[match(common, bz$record_id)]
    m <- if (length(common) && sum(av)>0 && sum(bv)>0)min(sum(av), sum(bv))else NA_real_
    for (stage in stages) {
      zz <- z[z$stage == stage,];
      eligible <- zz[is.finite(zz$count),]
      x <- eligible$count;
      c <- pd_concentration(x)
      qs <- if (length(x))unname(quantile(x, c(.25, .5, .75), type = 7))else rep(NA_real_, 3)
      positive <- is.finite(c['total']) && c['total']>0
      largest <- if (positive)paste(sort(eligible$record_id[x == max(x)]), collapse = ';')else NA_character_
      common_values <- if (stage == stages[1])av else bv
      si <- si + 1L;
      summary[[si]] <- data.frame(cohort, forest, species = sp, stage,
        source_records = length(ids), eligible_records = length(x), missing_count_records = nrow(zz) - length(x),
        total_listed_count = unname(c['total']), mean_per_record = if (length(x))mean(x)else NA_real_,
        median_per_record = qs[2], q1 = qs[1], q3 = qs[3], iqr = qs[3] - qs[1],
        minimum = if (length(x))min(x)else NA_real_, maximum = if (length(x))max(x)else NA_real_,
        records_with_listed_count = sum(x>0), record_frequency = if (length(x))mean(x>0)else NA_real_,
        record_frequency_percent = if (length(x))100 * mean(x>0)else NA_real_,
        mean_where_listed = if (any(x>0))mean(x[x>0])else NA_real_,
        largest_count_record_ids = largest, largest_record_share = unname(c['largest_record_share']),
        concentration = unname(c['concentration']),
        effective_number_of_survey_records = unname(c['effective_number_of_survey_records']),
        distinct_pair_concentration = unname(c['distinct_pair_concentration']),
        equal_count_common_records = length(common), equal_count_stage_total = sum(common_values),
        equal_count_target = m,
        expected_concentration_at_equal_count = pd_equal_count_expectation(common_values, m),
        location_independence_verified = FALSE)
      if (nrow(eligible)) {
        eligible <- eligible[order( - eligible$count, eligible$record_id),]
        ci <- ci + 1L;
        contributions[[ci]] <- data.frame(cohort, forest, species = sp, stage,
          record_id = eligible$record_id, listed_count = eligible$count,
          count_rank = rank( - eligible$count, ties.method = 'min'),
          count_share = if (positive)eligible$count / unname(c['total'])else NA_real_,
          cumulative_share = if (positive)cumsum(eligible$count) / unname(c['total'])else NA_real_,
          ordering = 'Descending count, record ID for ties; not geographical distance')
      }
    }
  }
  summary <- dplyr::bind_rows(summary);
  contributions <- dplyr::bind_rows(contributions)
  put(summary, 'abundance_frequency_and_concentration')
  put(summary[summary$cohort == 'primary',], 'primary_abundance_frequency_and_concentration')
  put(contributions, 'record_count_contributions')
  location <- q[, c('plot_id', 'forest', 'cohort', 'cultivated', 'primary_eligible',
      'lat', 'lon', 'coordinate_source', 'coordinate_verified', 'datum',
      'positional_uncertainty_m', 'boundary_verified', 'verified_area_m2')]
  location$mapping_eligible <- is.finite(location$lat) & is.finite(location$lon) &
  abs(location$lat) <= 90 & abs(location$lon) <= 180
  location$uncertainty_statement <- 'Worksheet coordinate; datum unconfirmed and accuracy unknown; not an individual plant position'
  put(location, 'worksheet_locations')
  derived <- file.path(root, 'analysis/Datasets/derived')
  dir.create(derived, recursive = TRUE, showWarnings = FALSE)
  write.csv(location, file.path(derived, 'worksheet_locations.csv'), row.names = FALSE, na = 'NA')
  put(location[!location$mapping_eligible,], 'mapping_exclusions')
  cols <- c(Ngangao = '#D4B58A', Mbololo = '#508A92')
  draw_counts <- function() {
    par(mfrow = c(2, 2), mar = c(6, 4.4, 3, 1))
    for (sp in species)for (stage in stages) {
      z <- ledger[ledger$primary_record & ledger$species == sp & ledger$stage == stage,]
      z <- z[order(match(z$forest, names(cols)), z$record_id),]
      finite <- is.finite(z$count);
      ymax <- max(c(1, z$count[finite])) * 1.15
      plot(seq_len(nrow(z)), z$count, type = 'h', lwd = 3, col = cols[z$forest], ylim = c(0, ymax),
        xaxt = 'n', xlab = '', ylab = 'Listed count per record',
        main = bquote(atop(italic(.(display(sp))), .(stage))), cex.main = .9)
      points(seq_len(nrow(z))[finite], z$count[finite], pch = 16, col = cols[z$forest[finite]], cex = .6)
      lab <- sub('Ngangao_', 'N', sub('Mbololo_', 'M', z$record_id))
      axis(1, at = seq_len(nrow(z)), labels = lab, las = 2, cex.axis = .6)
      mtext('Source record (N = Ngangao; M = Mbololo)', side = 1, line = 4, cex = .65)
      legend('topright', names(cols), pch = 16, col = cols, bty = 'n', cex = .7)
    }
  }
  draw_frequency <- function() {
    par(mfrow = c(1, 2), mar = c(7, 4.4, 3, 1))
    for (sp in species) {
      z <- summary[summary$cohort == 'primary' & summary$species == sp,]
      lab <- paste(z$forest, ifelse(z$stage == stages[1], 'Larger stems', 'Juveniles'), sep = '\n')
      mid <- barplot(z$record_frequency_percent, names.arg = lab, las = 2, cex.names = .8,
        col = cols[z$forest], ylim = c(0, 115), ylab = 'Recorded occurrence frequency (%)', main = bquote(italic(.(display(sp)))))
      text(mid, z$record_frequency_percent + 5, paste0(z$records_with_listed_count, '/', z$eligible_records), cex = .8)
    }
  }
  draw_concentration <- function() {
    par(mfrow = c(2, 2), mar = c(4.5, 4.4, 3, 1))
    for (sp in species)for (forest in names(cols)) {
      z <- contributions[contributions$cohort == 'primary' & contributions$species == sp & contributions$forest == forest,]
      n <- max(c(1, table(z$stage)))
      plot(0, 0, type = 'n', xlim = c(0, n), ylim = c(0, 1), xlab = 'Records ranked by count within each stage',
        ylab = 'Cumulative share of listed count', main = bquote(atop(italic(.(display(sp))), .(forest))), cex.main = .9)
      for (i in seq_along(stages)) {
        v <- z[z$stage == stages[i],]
        if (any(is.finite(v$cumulative_share)))lines(c(0, seq_len(nrow(v))), c(0, v$cumulative_share),
          type = 's', lwd = 2, lty = i, col = c('#0072B2', '#D55E00')[i])
      }
      legend('bottomright', vapply(stages, function(s) {
            v <- z[z$stage == s,];
            paste0(s, ' (N=', sum(v$listed_count), ')')
          }, character(1)), lty = 1:2, col = c('#0072B2', '#D55E00'), bty = 'n', cex = .65)
    }
  }
  draw_locator <- function() {
    par(mar = c(5, 5, 5, 1))
    z <- location[location$mapping_eligible,]
    if (!nrow(z)) {plot.new();
      text(.5, .5, 'No usable worksheet coordinates');
      return(invisible())}
    shape <- ifelse(z$cultivated, 15, ifelse(z$primary_eligible, 16, 17))
    plot(z$lon, z$lat, pch = shape, col = cols[z$forest], cex = 1.1,
      asp = 1 / cos(mean(z$lat) * pi / 180), xlab = 'Longitude (degrees)',
      ylab = 'Latitude (degrees)', main = 'Worksheet locations')
    mtext('Datum unconfirmed; positional accuracy unknown; symbols represent source records', side = 3, line = .5, cex = .75)
    legend('topright', c('Ngangao', 'Mbololo', 'Primary forest-associated', 'Identity unresolved', 'Cultivated'),
      pch = c(16, 16, 16, 17, 15), col = c(cols, 'black', 'black', 'black'), bty = 'n', cex = .8)
  }
  figures <- list(recorded_focal_abundance = draw_counts, record_frequency = draw_frequency,
    concentration_among_records = draw_concentration, provisional_record_locator = draw_locator)
  figure_names <- c(recorded_focal_abundance = 'Figure_3_recorded_focal_abundance',
    record_frequency = 'Figure_S5_recorded_occurrence_frequency',
    concentration_among_records = 'Figure_S6_concentration_among_survey_records',
    provisional_record_locator = 'Figure_S7_worksheet_source_record_locations')
  for (name in names(figures))for (ext in c('pdf', 'png')) {
    pd_temp_ready()
    width <- if (name == 'provisional_record_locator')9 else 12
    height <- if (name == 'record_frequency')6 else 9
    if (ext == 'pdf')pdf(rq_file(out, paste0(figure_names[[name]], '.pdf')), width = width, height = height)
    else png(rq_file(out, paste0(figure_names[[name]], '.png')), width = width * 250, height = height * 250, res = 250)
    tryCatch(figures[[name]](), finally = dev.off())
  }
  put(data.frame(method = c('Density extrapolation', 'Detection-corrected occupancy',
        'Ripley K and pair correlation', 'Habitat-loss redundancy', 'Demographic trends'), status = 'Deferred',
      reason = c('Unknown footprints, overlap and regional inclusion probabilities',
        'No repeat detection histories', 'No verified individual coordinates or observation window',
        'Record concentration is not independent threatened-population representation',
        'No linked repeated census')), 'deferred_methods')
  notes <- c('# Recorded abundance and distribution methods', '',
    'Objective: describe recorded abundance, occurrence frequency and concentration among survey records of Coffea fadenii and Afrocarpus usambarensis, separately by forest and operational stage.',
    'This analysis reads the existing vegetation-structure eligibility ledger, stem audit and taxonomy decisions, together with the original-workbook juvenile row audit generated by 01_woody_composition.R. No missing configuration or reporting-helper folders are required. Primary forest-associated records exclude two cultivated Mbololo cases and three identity-flagged records. Broader forest-context and cultivated cases remain separate descriptive cohorts. All 28 source records are retained in the ledger. Accepted-name decisions are copied from the reconciled input tables; unresolved labels are not silently merged with focal taxa.',
    'Counts are larger-stem source rows or sums of listed juvenile-table counts. Headers remain a separate evidence stream; no hypothetical header addition enters any listed-count summary. An absent matching focal row gives a listed zero only when stage-table evidence exists. When no stage-table rows were extracted, the count is NA, and its record is excluded from that stage denominator. Recorded zeros do not certify ecological absence or complete detection.',
    'For every cohort, species, forest and stage, report source and eligible record numbers, missing-count records, total, mean, median, type-7 quartiles, IQR and range. Mean where listed is conditional on positive counts and is undefined when no record has a positive count. Totals describe listed observations, not unique genets or whole-forest populations.',
    'Record frequency is the number of eligible records with a positive listed count divided by eligible records examined. Both numerator and denominator are reported. This is not detection-corrected occupancy. No binomial uncertainty interval or significance test assumes independent random records.',
    'For total N>0, count shares p_i=n_i/N give empirical concentration C=sum(p_i^2), effective number of survey records 1/C and largest-record share max(p_i). Largest-record IDs and ties are retained. N=0 makes shares, concentration, effective number and largest share undefined, rather than zero. An entirely missing stratum has total NA.',
    'Distinct-pair concentration is sum(n_i*(n_i-1))/(N*(N-1)), defined for N>=2. Stage comparisons additionally use only records with available counts for both stages and m=min(N_larger,N_juvenile)>0. Expected empirical concentration after drawing m observed individuals without replacement is 1/m+(m-1)*C_pair/m. For m=1 it equals one. This exact arithmetic subsampling expectation is not an ecological null or a correction for unequal effort. If either common-record stage total is zero, the equal-count comparison is undefined. Common record numbers and stage totals are supplied beside the calculation.',
    'Cumulative count shares rank records separately for each stage; lines do not connect geographical neighbours or demonstrate demographic replacement. Latitude/longitude symbols identify source records and context, not individual plants. Datum and positional accuracy remain unconfirmed; uncertainty circles are not invented. Missing or invalid coordinates are excluded from the locator and retained in a mapping-exclusion table. No forest boundaries, density surfaces or habitat suitability are inferred.',
    'No population-density extrapolation, geographical clustering test, survival estimate, recovery assessment or redundancy against habitat loss is performed. Cultivated records cannot establish wild occurrence. Aggregation statistics depend on source inclusion and observation processes.',
    '', 'Published precedents and explicit adaptations',
    'Buckley and Freckleton (2010), Journal of Ecology 98, 645-658: distinguish local plant incidence from occupancy across sampling units. https://doi.org/10.1111/j.1365-2745.2010.01650.x',
    'Pinho et al. (2022), Diversity and Distributions 28, 1497-1511: distinguish local abundance from occurrence across explicitly defined units and scales. Our record-based descriptive summaries do not reproduce their regional models. https://doi.org/10.1111/ddi.13587',
    'Perret et al. (2023), Journal of Ecology 111, 1483-1496: imperfect and heterogeneous detection limits interpretation of raw plant counts. Detection is not estimated from these one-off records. https://doi.org/10.1111/1365-2745.14110',
    'The Simpson-type mathematics is adapted to count shares among records; it must not be presented as species diversity or individual point-pattern analysis. Count-standardised expectations are algebraic descriptions of the observed counts.',
    '', 'Run: Rscript --vanilla "analysis/scripts/03_focal_abundance.R"',
    'Outputs: analysis/outputs/research question 3 (figures/ and tables/). Existing manuscripts and analysis scripts are not edited.')
  writeLines(notes, rq_file(out, 'METHODS.md'))
  results <- c('# Research question 3: Recorded abundance, frequency and concentration', '',
    'Primary summaries concern selected source records; sampling area, individual identity and detection remain unverified.')
  primary <- summary[summary$cohort == 'primary',]
  metrics <- c('total_listed_count', 'median_per_record', 'record_frequency_percent', 'concentration')
  compare <- function(group, high, low, by) {
    a <- primary[primary[[group]] == high, c(by, 'eligible_records', metrics)]
    b <- primary[primary[[group]] == low, c(by, 'eligible_records', metrics)]
    d <- merge(a, b, by = by, suffixes = c('_high', '_low'))
    d$high_group <- high;
    d$low_group <- low
    for (m in metrics)d[[paste0(m, '_difference')]] <- d[[paste0(m, '_high')]] - d[[paste0(m, '_low')]]
    d
  }
  put(compare('forest', 'Mbololo', 'Ngangao', c('species', 'stage')), 'primary_forest_contrasts')
  put(compare('stage', 'Recorded juveniles', 'Listed larger stems', c('species', 'forest')), 'primary_size_class_contrasts')
  results <- c(results, 'Descriptive forest and size-class differences are tabulated alongside their eligible-record denominators. Raw count and concentration differences remain sensitive to effort and count totals.')
  for (i in seq_len(nrow(primary))) {
    z <- primary[i,]
    results <- c(results, sprintf('%s; %s; %s: total %.0f, median %.2f [Q1 %.2f, Q3 %.2f]; recorded in %d/%d eligible records (%.1f%%).',
        display(z$species), z$forest, z$stage, z$total_listed_count, z$median_per_record, z$q1, z$q3,
        z$records_with_listed_count, z$eligible_records, z$record_frequency_percent))
    if (is.finite(z$total_listed_count) && z$total_listed_count>0)results <- c(results, sprintf('Largest-record share %.1f%% (%s); concentration %.4f; effective number of survey records %.2f.',
        100 * z$largest_record_share, z$largest_count_record_ids, z$concentration, z$effective_number_of_survey_records))
  }
  results <- c(results, '', 'No listed count is not verified absence. Concentration is among records, not geographical clustering or habitat-loss redundancy. Cultivated plants, missing stage evidence and focal headers are reported separately.')
  writeLines(results, rq_file(out, 'RESULTS.md'))
  population_reporting(out, primary, summary, ledger, headers, location)
  inputs <- c(file.path(here, '03_focal_abundance.R'), source_paths)
  put(data.frame(path = substring(inputs, nchar(root) + 2), md5 = unname(tools::md5sum(inputs))),
    'input_and_code_hashes')
  writeLines(capture.output(sessionInfo()), rq_file(out, 'sessionInfo.txt'))
  jsonlite::write_json(list(primary_records = sum(q$primary_eligible), source_records = nrow(q),
      species = species, stages = stages, density_estimated = FALSE, occupancy_detection_corrected = FALSE,
      spatial_process_tests = FALSE, coordinates_verified = FALSE, confirmatory_tests = FALSE,
      count_standardisation = 'Exact expected empirical concentration under equal-total subsampling'),
    rq_file(out, 'analysis_status.json'), auto_unbox = TRUE, pretty = TRUE)
  message('Recorded population distribution complete: ', out)
  two_part <- pd_two_part_report(ledger, ledger$primary_record, out)
  invisible(list(ledger = ledger, summary = summary, two_part = two_part, contributions = contributions, output = out))
}

# ============================================================
# 9. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep('^--file=', commandArgs(), value = TRUE)
  start <- if (length(arg))dirname(normalizePath(sub('^--file=', '', arg[1]), mustWork = TRUE))else getwd()
  root <- pd_project_root(start)
  run_population_distribution(root, file.path(root, 'analysis/scripts'))
  source(file.path(file.path(root, 'analysis'), 'scripts', 'combine_output_tables.R'), local = TRUE)
  combine_output_tables(file.path(file.path(root, 'analysis'), 'outputs', 'research question 3'), '03_focal_abundance.R')
}
