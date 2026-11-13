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
# 1. Floristic temp ready
# ============================================================

# Research question 1: How does recorded woody taxon composition differ between the surveyed neighbourhoods of Ngangao and Mbololo when juveniles and larger stems are considered separately?
# Floristic composition of surveyed woody neighbourhoods in the Taita Hills.
# density scaling, rarefaction or extrapolation is performed.
floristic_temp_ready <- function() {
  tmp <- tempdir(check = TRUE)
  if (!dir.exists(tmp))dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  probe <- tempfile('floristic_write_check_', tmpdir = tmp)
  writable <- suppressWarnings(file.create(probe))
  if (isTRUE(writable))unlink(probe)
  if (!dir.exists(tmp) || !isTRUE(writable))
  stop('R temporary directory is unavailable or not writable: ', tmp,
    '\nRestart the R session and Source this script again.', call. = FALSE)
  invisible(tmp)
}

# ============================================================
# 2. Floristic project root
# ============================================================

floristic_project_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = '/', mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, 'analysis/scripts/01_woody_composition.R'))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  stop('Cannot locate the project root from: ', start,
    '\nSet the working directory to Franklinia.Taita Taveta.Project or its analysis/scripts directory.', call. = FALSE)
}

# ============================================================
# 3. Floristic provisional
# ============================================================

floristic_provisional <- function(root, out, export = TRUE) {
  rows <- list();
  records <- list();
  evidence <- list()
  norm <- function(x) tolower(trimws(gsub('[[:space:]\u00a0]+', ' ', x, perl = TRUE)))
  for (forest in c('Ngangao', 'Mbololo')) {
    rel <- file.path('analysis/Datasets/field', paste0(forest, '.Mapping.Data.Franklinia.xlsx'))
    path <- file.path(root, rel)
    for (sh in readxl::excel_sheets(path)) {
      id <- as.integer(sub('Plot.', '', sh, fixed = TRUE))
      stopifnot(!is.na(id))
      key <- sprintf('%s_%02d', forest, id)
      reason <- if (forest != 'Mbololo') '' else if (id %in% 1:2)
      'Roadside/farm context; conservatively excluded from forest comparison' else if (id %in% c(3, 5, 6))
      'Identity exclusion documented in analysis/Datasets/field/Manuscript_strengthening_plan.md' else ''
      records[[key]] <- data.frame(record_id = key, forest, workbook = rel, sheet = sh, included = !nzchar(reason), reason)
      d <- readxl::read_excel(path, sheet = sh, range = readxl::cell_limits(c(1, 1), c(NA, 8)),
        col_names = FALSE, col_types = 'text', .name_repair = 'minimal')
      field <- function(pattern) {
        cells <- unlist(d[seq_len(min(10L, nrow(d))),, drop = FALSE], use.names = FALSE)
        hits <- cells[which(grepl(pattern, trimws(cells), ignore.case = TRUE))]
        if (length(hits) == 1L)trimws(sub(pattern, '', trimws(hits), ignore.case = TRUE))else NA_character_
      }
      records[[key]]$worksheet_plot_number <- id
      records[[key]]$header_plot_number <- suppressWarnings(as.integer(field('^plot:[[:space:]]*')))
      records[[key]]$header_forest <- field('^forest:[[:space:]]*')
      records[[key]]$survey_date_label <- field('^date:[[:space:]]*')
      records[[key]]$focal_header_taxon <- norm(field('^name of the species:[[:space:]]*'))
      records[[key]]$focal_header_dbh <- field('^dbh:[[:space:]]*')
      records[[key]]$focal_header_saplings <- field('^no of saplings:[[:space:]]*')
      records[[key]]$focal_header_larger_count <- field('^no of mature( trees)?:[[:space:]]*')
      records[[key]]$latitude_label <- field('^latitude:[[:space:]]*')
      records[[key]]$longitude_label <- field('^longitude:[[:space:]]*')
      records[[key]]$worksheet_context <- field('^description of the plot and vegetation:[[:space:]]*')
      for (stage in c('Listed larger stems', 'Recorded juveniles')) {
        juvenile <- stage == 'Recorded juveniles'
        cols <- if (juvenile)c(6, 7, 8)else c(1, 2, 3)
        h <- which(norm(d[[cols[2]]]) == 'species name' &
          grepl(if (juvenile)'number' else 'dbh', norm(d[[cols[3]]])))
        e <- paste(key, stage)
        evidence[[e]] <- data.frame(record_id = key, stage, status = 'Survey status unknown',
          header_row = if (length(h) == 1L)h else NA_integer_, numbered_rows = 0L,
          usable_rows = 0L, invalid_rows = 0L, listed_total = NA_real_,
          count_or_dbh_header = if (length(h) == 1L)d[[cols[3]]][h]else NA_character_,
          height_header = if (length(h) == 1L && !juvenile)d[[4]][h]else NA_character_)
        if (length(h)>1L)stop('Ambiguous ', stage, ' header in ', key)
        if (!length(h))next
        ix <- which(seq_len(nrow(d))>h)
        stop_at <- which(norm(d[[cols[1]]][ix]) == 'no' |
          grepl('^(shrubs?|herbs?|herbaceous)', norm(d[[cols[2]]][ix])))
        if (length(stop_at))ix <- head(ix, min(stop_at) - 1L)
        ix <- ix[which(grepl('^[0-9]+$', trimws(d[[cols[1]]][ix])) &
            !is.na(d[[cols[2]]][ix]) & nzchar(trimws(d[[cols[2]]][ix])))]
        count <- if (juvenile)suppressWarnings(as.numeric(d[[8]][ix]))else rep(1, length(ix))
        usable <- is.finite(count) & count >= 0 & count == floor(count)
        evidence[[e]]$numbered_rows <- length(ix)
        evidence[[e]]$usable_rows <- sum(usable)
        evidence[[e]]$invalid_rows <- sum(!usable)
        if (length(ix) && all(usable))evidence[[e]]$listed_total <- sum(count)
        if (length(ix))evidence[[e]]$status <- if (any(!usable))'Invalid counts' else
        if (sum(count)>0)'Recorded counts available'else 'Surveyed: explicit zero counts'
        if (length(ix)) rows[[paste(key, stage)]] <- data.frame(record_id = key, forest, stage,
          workbook = rel, sheet = sh, source_row = ix, source_label = d[[cols[2]]][ix],
          recorded_label = norm(d[[cols[2]]][ix]), count,
          raw_dbh = if (juvenile)NA_character_ else d[[3]][ix],
          raw_height = if (juvenile)NA_character_ else d[[4]][ix],
          usable_count = is.finite(count) & count >= 0 & count == floor(count))
      }
    }
  }
  rec <- do.call(rbind, records);
  raw <- do.call(rbind, rows)
  evidence <- do.call(rbind, evidence)
  same_coord <- !is.na(rec$latitude_label) & !is.na(rec$longitude_label)
  key <- paste(rec$latitude_label, rec$longitude_label, sep = '|')
  rec$coordinate_label_shared <- same_coord & (duplicated(key) | duplicated(key, fromLast = TRUE))
  rec$header_plot_matches_sheet <- rec$header_plot_number == rec$worksheet_plot_number
  rec$header_forest_matches_workbook <- tolower(rec$header_forest) == tolower(rec$forest)
  rec$footprint_verified <- FALSE
  rec$coordinate_datum_verified <- FALSE
  rec$juvenile_search_area_verified <- FALSE
  rec$larger_stem_search_area_verified <- FALSE
  rec$record_independence_verified <- FALSE
  rownames(rec) <- rownames(raw) <- rownames(evidence) <- NULL
  if (!export)return(list(records = rec, source_rows = raw, evidence = evidence))
  put <- function(x, n)write.csv(x, rq_file(out, paste0(n, '.csv')), row.names = FALSE)
  put(rec, 'worksheet_selection');
  put(raw, 'source_rows');
  put(evidence, 'worksheet_stage_evidence_audit')
  floristic_composition_report(rec, raw, evidence, out, provisional = TRUE, root = root)
}

# ============================================================
# 4. Floristic inputs
# ============================================================

floristic_inputs <- function(root, out) {
  paths <- rq_input(root, c('analysis/population_distribution/record_reconciliation.csv',
      'analysis/vegetation_structure/taxonomic_decisions.csv',
      'analysis/vegetation_structure/individual_stem_measurement_audit.csv'))
  missing <- !file.exists(paths)
  if (any(missing)) {
    status <- c('Audited Question 1 unavailable: required survey audits are missing.',
      'Restore the original audit CSVs to these locations:',
      paste0(' - ', paths[missing]),
      'Place independently reconciled inputs in analysis/Datasets/audits.',
      'The field workbooks contain raw observations; they do not replace the audited record eligibility, woody taxon decisions, and stem-row exclusions.',
      'Running a separate provisional worksheet comparison; no audited publication tables are regenerated.')
    message(paste(status, collapse = '\n'))
    return(NULL)
  }
  q <- read.csv(paths[1], stringsAsFactors = FALSE, check.names = FALSE)
  tax <- read.csv(paths[2], stringsAsFactors = FALSE, check.names = FALSE)
  stems <- read.csv(paths[3], stringsAsFactors = FALSE, check.names = FALSE)
  norm <- function(x)tolower(trimws(gsub('[[:space:]\u00a0]+', ' ', as.character(x), perl = TRUE)))
  tax$source_name <- norm(tax$source_name);
  tax$analysis_name <- norm(tax$analysis_name)
  stopifnot(!anyDuplicated(tax$source_name))
  taxa <- sort(unique(tax$analysis_name[tax$woody_eligible]))
  A <- J <- matrix(0, nrow(q), length(taxa), dimnames = list(q$plot_id, taxa))
  for (i in seq_len(nrow(stems)))if (isTRUE(stems$woody_eligible[i])) {
    row <- match(stems$plot_id[i], rownames(A));
    col <- match(stems$analysis_name[i], colnames(A))
    stopifnot(!is.na(row), !is.na(col));
    A[row, col] <- A[row, col] + 1
  }
  juvenile_rows <- list();
  k <- 0L;
  workbooks <- character()
  for (i in seq_len(nrow(q))) {
    # Resolve portable paths using retained workbook names, not another machine's absolute root.
    workbook <- file.path(root, 'analysis/Datasets/field', basename(gsub('\\', '/', q$source_file[i], fixed = TRUE)))
    if (!file.exists(workbook))stop('Original workbook missing: ', workbook, call. = FALSE)
    workbooks <- c(workbooks, workbook)
    dat <- suppressMessages(readxl::read_excel(workbook, sheet = q$sheet[i], col_names = FALSE, .name_repair = 'minimal'))
    if (ncol(dat)<8)stop('Juvenile source table cannot be located: ', q$plot_id[i])
    headers <- which(grepl('species name', tolower(as.character(dat[[7]]))) & grepl('number', tolower(as.character(dat[[8]]))))
    if (length(headers) != 1)stop('Expected one juvenile-table header: ', q$plot_id[i])
    # A separate shrub/herb list below the juvenile table is not a count table.
    stops <- which(seq_len(nrow(dat))>headers & (tolower(trimws(as.character(dat[[6]]))) == 'no' | grepl('^(shrubs?|herbs?|herbaceous)', norm(dat[[7]]))))
    last <- if (length(stops))min(stops) - 1L else nrow(dat)
    for (rr in seq.int(headers + 1L, last)) {
      label <- norm(dat[[7]][rr]);
    value <- as.character(dat[[8]][rr])
      if (is.na(label) || !nzchar(label))next
      count <- suppressWarnings(as.numeric(value));
      m <- match(label, tax$source_name)
      eligible <- !is.na(m) && isTRUE(tax$woody_eligible[m])
      if (eligible && (!is.finite(count) || count<0 || count != floor(count)))
      stop('Unresolved juvenile count in ', q$plot_id[i], ', worksheet row ', rr, ': ', value)
      k <- k + 1L;
      juvenile_rows[[k]] <- data.frame(plot_id = q$plot_id[i], source_file = workbook,
        sheet = q$sheet[i], source_row = rr, source_name = label, analysis_name = if (!is.na(m))tax$analysis_name[m]else label,
        count = count, woody_eligible = eligible, decision = if (is.na(m))'Unresolved label; excluded from woody matrix'else tax$decision[m])
      if (eligible)J[i, match(tax$analysis_name[m], colnames(J))] <- J[i, match(tax$analysis_name[m], colnames(J))] + count
    }
  }
  juvenile_audit <- dplyr::bind_rows(juvenile_rows)
  write.csv(juvenile_audit, rq_file(out, 'juvenile_source_row_audit.csv'), row.names = FALSE)
  write.csv(juvenile_audit, rq_input(root, 'juvenile_source_row_audit.csv'), row.names = FALSE)
  write.csv(data.frame(record_id = rownames(A), A, check.names = FALSE), rq_file(out, 'larger_stem_matrix.csv'), row.names = FALSE)
  write.csv(data.frame(record_id = rownames(J), J, check.names = FALSE), rq_file(out, 'juvenile_matrix.csv'), row.names = FALSE)
  # Independent focal ledger and stem totals detect extraction or mapping changes.
  check <- q$primary_eligible
  stopifnot(sum(A[check,, drop = FALSE]) == sum(stems$woody_eligible & stems$plot_id %in% q$plot_id[check]),
    all(J[, 'coffea fadenii'] == q$juvenile_cf_listed))
  list(q = q, tax = tax, A = A, J = J, files = c(paths, unique(workbooks)))
}

# ============================================================
# 5. Floristic reference
# ============================================================

floristic_reference <- function(root) {
  h <- readxl::read_excel(file.path(root, 'analysis/Datasets/field/Data.Taita.xlsx'),
    sheet = 'Species.list', skip = 2, col_names = FALSE, col_types = 'text', .name_repair = 'minimal')
  norm <- function(x)tolower(trimws(gsub('[[:space:]\u00a0]+', ' ', x, perl = TRUE)))
  data.frame(taxon = norm(h[[3]]), family = trimws(h[[2]]), habit = toupper(trimws(h[[4]])))
}

# ============================================================
# 6. Floristic wcvp families
# ============================================================

# Match the taxa retained in this run directly against the installed WCVP snapshot.
floristic_wcvp_families <- function(recorded, workbook_families, out, derived) {
  queries <- data.frame(taxon = recorded, query = paste0(toupper(substr(recorded, 1, 1)), substring(recorded, 2)))
  checklist <- rWCVPdata::wcvp_names
  matches <- rWCVP::wcvp_match_exact(queries, checklist, name_col = 'query', id_col = 'taxon')
  accepted_id <- ifelse(matches$wcvp_status %in% 'Accepted', matches$wcvp_id, matches$wcvp_accepted_id)
  at <- match(accepted_id, checklist$plant_name_id)
  valid <- !is.na(at) & checklist$taxon_status[at] %in% 'Accepted'
  matches$accepted_id <- ifelse(valid, accepted_id, NA_real_)
  matches$accepted_name <- ifelse(valid, checklist$taxon_name[at], NA_character_)
  matches$accepted_authors <- ifelse(valid, checklist$taxon_authors[at], NA_character_)
  matches$wcvp_family <- ifelse(valid, checklist$family[at], NA_character_)
  matches$powo_id <- ifelse(valid, checklist$powo_id[at], NA_character_)
  # Distribution evidence narrows the nomenclatural candidates without confirming a field identification.
  distribution <- rWCVPdata::wcvp_distributions
  kenya_native <- unique(distribution$plant_name_id[which(distribution$area_code_l3 == 'KEN' &
        distribution$introduced == 0 & distribution$extinct == 0 & distribution$location_doubtful == 0)])
  matches$native_kenya_in_wcvp <- ifelse(is.na(matches$accepted_id), NA, matches$accepted_id %in% kenya_native)
  resolved <- do.call(rbind, lapply(split(matches, matches$taxon), function(m) {
        t <- m$taxon[1]
        hits <- unique(m$accepted_id[!is.na(m$accepted_id)])
        chosen <- if (length(hits) == 1L)match(hits, m$accepted_id)else NA_integer_
        families <- unique(m$wcvp_family[!is.na(m$wcvp_family) & nzchar(m$wcvp_family)])
        f <- if (length(families) == 1L)families else NA_character_
        ok <- !is.na(f) && nzchar(f)
        local <- workbook_families[[t]]
        status <- if (ok && length(hits) == 1L) 'Resolved' else if (ok) 'Family resolved; accepted taxon ambiguous' else if (all(is.na(m$wcvp_id))) 'No exact match' else
        if (length(hits)>1L)'Ambiguous accepted taxon'else 'Accepted family unresolved'
        data.frame(recorded_taxon = t, accepted_name = m$accepted_name[chosen],
          accepted_authors = m$accepted_authors[chosen], wcvp_family = f,
          workbook_family = paste(local, collapse = ' | '), resolution = status,
          family_comparison = if (ok)if (f %in% local)'Same'else 'Different'else NA_character_,
          exact_matches = sum(!is.na(m$wcvp_id)),
          kenya_supported_candidate = if (length(hits)>1L) {
            n <- unique(m$accepted_name[which(m$native_kenya_in_wcvp)])
            if (length(n) == 1L)n else NA_character_
          } else NA_character_,
          accepted_wcvp_id = m$accepted_id[chosen],
          powo_id = m$powo_id[chosen], rWCVP_version = as.character(packageVersion('rWCVP')),
          rWCVPdata_version = as.character(packageVersion('rWCVPdata')))
      }))
  stopifnot(nrow(resolved) == length(recorded), !anyDuplicated(resolved$recorded_taxon))
  write.csv(resolved, file.path(derived, 'WCVP_species_family_lookup.csv'), row.names = FALSE, na = '')
  write.csv(resolved, rq_file(out, 'WCVP_species_family_lookup.csv'), row.names = FALSE, na = '')
  write.csv(matches, file.path(derived, 'WCVP_exact_match_audit.csv'), row.names = FALSE, na = '')
  write.csv(matches, rq_file(out, 'WCVP_exact_match_audit.csv'), row.names = FALSE, na = '')
  resolved
}

# ============================================================
# 7. Floristic composition report
# ============================================================

floristic_composition_report <- function(rec, rows, evidence, out, provisional, repeats = 8L, matrices = NULL, root = NULL) {
  prefix <- if (provisional)'provisional_'else ''
  put <- function(x, n, fileEncoding = '')write.csv(x, rq_file(out, paste0(n, '.csv')),
    row.names = FALSE, na = 'NA', fileEncoding = fileEncoding)
  forests <- c('Ngangao', 'Mbololo')
  stages <- c('Recorded juveniles', 'Listed larger stems')
  cols <- setNames(c('#D4B58A', '#508A92'), forests)
  symbols <- setNames(c(16, 17), forests)
  if (provisional) {
    habit_source <- 'analysis/Datasets/field/Data.Taita.xlsx'
    h <- readxl::read_excel(file.path(root, habit_source), sheet = 'Species.list', col_names = FALSE,
      col_types = 'text', .name_repair = 'minimal')
    norm <- function(x)tolower(trimws(gsub('[[:space:]\u00a0]+', ' ', x, perl = TRUE)))
    names_h <- norm(h[[3]])
    codes <- toupper(trimws(h[[4]]))
    labels <- sort(unique(rows$recorded_label))
    decision <- do.call(rbind, lapply(labels, function(x) {
          c <- sort(unique(codes[which(names_h == x & !is.na(codes) & nzchar(codes))]))
          include <- length(c)>0 && all(grepl('^[TSL](/[TSL])*$', c))
          data.frame(recorded_label = x, source_habit = paste(c, collapse = '; '), included_as_woody = include,
            decision = if (include)'Documented tree, shrub or liana'else if (length(c))
            'Habit not unequivocally woody'else 'No exact habit match', source = habit_source)
        }))
    put(decision, 'woody_habit_screen')
    rows <- rows[rows$recorded_label %in% decision$recorded_label[decision$included_as_woody],]
  }
  selected <- rec[rec$included,]
  e <- merge(evidence, rec[, c('record_id', 'forest', 'included')], by = 'record_id')
  reference <- floristic_reference(root)
  derived <- file.path(root, 'analysis/Datasets/derived')
  dir.create(derived, recursive = TRUE, showWarnings = FALSE)
  stage_inputs <- lapply(seq_along(stages), function(k) {
      stage <- stages[k]
      info <- e[e$stage == stage & e$included,]
      info <- info[match(selected$record_id, info$record_id),]
      stopifnot(identical(info$record_id, selected$record_id))
      raw <- rows[rows$stage == stage & rows$record_id %in% info$record_id & rows$usable_count,]
      if (is.null(matrices)) {
        labels <- sort(unique(raw$recorded_label))
        counts <- matrix(0, nrow(info), length(labels), dimnames = list(info$record_id, labels))
        for (i in seq_len(nrow(raw)))counts[raw$record_id[i], raw$recorded_label[i]] <-
        counts[raw$record_id[i], raw$recorded_label[i]] + raw$count[i]
      } else counts <- matrices[[k]][info$record_id,, drop = FALSE]
      observed <- info$status %in% c('Recorded counts available', 'Surveyed: explicit zero counts')
      info$status[observed & rowSums(counts) == 0 & info$status == 'Recorded counts available'] <- 'Surveyed: no retained counts'
      counts[!observed,] <- NA_real_
      list(info = info, counts = counts, observed = observed)
    })
  recorded <- sort(unique(unlist(lapply(stage_inputs, function(z)
          colnames(z$counts)[colSums(z$counts, na.rm = TRUE)>0]))))
  stopifnot(length(recorded)>0)
  workbook_families <- vapply(recorded, function(t) {
      f <- unique(reference$family[which(reference$taxon == t)])
      f <- f[!is.na(f) & nzchar(f)]
      if (length(f) == 1L)f else 'Unresolved family'
    }, character(1))
  lookup <- floristic_wcvp_families(recorded, workbook_families, out, derived)
  families <- richness <- occurrences <- family_audit <- list()
  summaries <- pairs <- scores <- diagnostics <- taxa <- availability <- list()
  for (k in seq_along(stages)) {
    stage <- stages[k]
    tag <- c('juveniles', 'larger_stems')[k]
    info <- stage_inputs[[k]]$info
    counts <- stage_inputs[[k]]$counts
    observed <- stage_inputs[[k]]$observed
    stage_counts <- data.frame(record_id = info$record_id, forest = info$forest, status = info$status, counts, check.names = FALSE)
    display_counts <- stage_counts
    names(display_counts) <- sub('^Ocotea usambarensis$', 'Kuloa usambarensis',
      names(stage_counts), ignore.case = TRUE)
    put(display_counts, paste0(tag, '_record_taxon_counts'))
    if (k == 1L)write.csv(stage_counts, file.path(derived, 'juveniles_record_taxon_counts.csv'), row.names = FALSE, na = 'NA')
    use <- observed & rowSums(counts, na.rm = TRUE)>0
    info$included_in_composition <- use
    put(info, paste0(tag, '_record_availability'))
    availability[[k]] <- do.call(rbind, lapply(forests, function(f) {
          z <- info[info$forest == f,]
          data.frame(stage, forest = f, selected_records = nrow(z), nonempty = sum(z$included_in_composition),
            surveyed_empty = sum(z$status == 'Surveyed: explicit zero counts'),
            surveyed_no_retained_woody_counts = sum(z$status == 'Surveyed: no retained counts'),
            never_surveyed = sum(z$status == 'Not surveyed'), survey_status_unknown = sum(z$status == 'Survey status unknown'),
            invalid_counts = sum(z$status == 'Invalid counts'))
        }))
    x <- counts[use,, drop = FALSE]
    group <- info$forest[use]
    if (nrow(x)<2L || !all(forests %in% group))stop('Insufficient nonempty records in both forests for ', stage)
    labels <- colnames(x)[colSums(x)>0]
    workbook_family <- unname(workbook_families[labels])
    at <- match(labels, lookup$recorded_taxon)
    family <- lookup$wcvp_family[at]
    family[is.na(family) | !nzchar(family)] <- 'Unresolved family'
    resolution <- lookup$resolution[at]
    resolution[is.na(at)] <- 'Recorded label absent from WCVP lookup'
    family_audit[[k]] <- data.frame(stage, taxon = labels, family, workbook_family,
      accepted_name = lookup$accepted_name[at], resolution,
      kenya_supported_candidate = if ('kenya_supported_candidate' %in% names(lookup))
      lookup$kenya_supported_candidate[at] else NA_character_,
      source = 'WCVP exact name matching; family of accepted taxon or consensus among candidates')
    totals <- sapply(forests, function(f)colSums(x[group == f, labels, drop = FALSE]))
    occurrences[[k]] <- data.frame(stage, taxon = labels, family,
      Ngangao_stems = totals[, 1], Mbololo_stems = totals[, 2],
      recorded_in = ifelse(rowSums(totals>0) == 2, 'Both forests',
        ifelse(totals[, 1]>0, 'Ngangao only', 'Mbololo only')))
    richness[[k]] <- data.frame(stage, forest = forests, nonempty_records = as.integer(table(factor(group, levels = forests))),
      observed_taxa = colSums(totals>0), shared_taxa = sum(rowSums(totals>0) == 2),
      forest_only_taxa = colSums(totals>0) - sum(rowSums(totals>0) == 2))
    families[[k]] <- do.call(rbind, lapply(forests, function(f) {
          n <- totals[, f];
          assigned <- family != 'Unresolved family'
          taxa_denom <- sum(n>0 & assigned);
          stem_denom <- sum(n[assigned])
          do.call(rbind, lapply(sort(unique(family)), function(g) {
                ix <- family == g
                data.frame(stage, forest = f, family = g, recorded_taxa = sum(n[ix]>0), recorded_stems = sum(n[ix]),
                  all_recorded_taxa = sum(n>0), all_recorded_stems = sum(n),
                  assigned_taxa_denominator = taxa_denom, assigned_stems_denominator = stem_denom,
                  taxon_share_percent = if (g != 'Unresolved family' && taxa_denom)100 * sum(n[ix]>0) / taxa_denom else NA_real_,
                  stem_share_percent = if (g != 'Unresolved family' && stem_denom)100 * sum(n[ix]) / stem_denom else NA_real_)
              }))
        }))
    rel <- sweep(x, 1, rowSums(x), '/')
    stopifnot(all(abs(rowSums(rel) - 1)<1e-10))
    means <- sapply(forests, function(f)colMeans(rel[group == f,, drop = FALSE]))
    td <- data.frame(stage, taxon = colnames(rel), Ngangao_percent = 100 * means[, 1], Mbololo_percent = 100 * means[, 2])
    td$difference_pp <- td$Mbololo_percent - td$Ngangao_percent
    td <- td[order( - abs(td$difference_pp), td$taxon),]
    td$displayed <- seq_len(nrow(td)) <= 10L
    taxa[[k]] <- td
    for (method in c('Bray-Curtis', 'Jaccard')) {
      d <- if (method == 'Bray-Curtis')vegan::vegdist(rel, method = 'bray') else vegan::vegdist(x, method = 'jaccard', binary = TRUE)
      stopifnot(all(is.finite(d)), all(d >= 0 & d <= 1 + 1e-10))
      ij <- which(lower.tri(as.matrix(d)), arr.ind = TRUE)
      typ <- ifelse(group[ij[, 1]] == group[ij[, 2]], paste('Within', group[ij[, 1]]), 'Between fragments')
      pair <- data.frame(stage, method, record_1 = rownames(x)[ij[, 1]], record_2 = rownames(x)[ij[, 2]], comparison = typ, dissimilarity = as.numeric(d))
      id <- paste(k, method)
      pairs[[id]] <- pair
      summaries[[id]] <- do.call(rbind, lapply(c('Between fragments', paste('Within', forests)), function(g) {
            v <- pair$dissimilarity[pair$comparison == g]
            q <- if (length(v))unname(quantile(v, c(.25, .5, .75)))else rep(NA_real_, 3)
            data.frame(stage, method, comparison = g, n_pairs = length(v), median = q[2], q1 = q[1], q3 = q[3], iqr = q[3] - q[1])
          }))
      if (nrow(x)<6L || length(unique(as.numeric(d)))<2L) {
        diagnostics[[id]] <- data.frame(stage, method, attempt = 0, stress = NA_real_, repeated_best = NA_integer_, selected = FALSE,
          message = 'Insufficient records or nondegenerate distances for 2D NMDS')
        next
      }
      fits <- vector('list', repeats);
      checks <- vector('list', repeats)
      for (i in seq_len(repeats)) {
        set.seed(5162026 + k * 100L + if (method == 'Jaccard')i + 20L else i)
        notes <- character()
        fit <- tryCatch(withCallingHandlers(vegan::metaMDS(d, k = 2, try = 20, trymax = 60, autotransform = FALSE,
              wascores = FALSE, trace = 0, maxit = 300), warning = function(w) {
              notes <<- c(notes, conditionMessage(w));
              invokeRestart('muffleWarning')
            }), error = function(err) {notes <<- c(notes, conditionMessage(err));
            NULL})
        fits[i] <- list(fit)
        checks[[i]] <- data.frame(stage, method, attempt = i, stress = if (is.null(fit))NA_real_ else fit$stress,
          repeated_best = if (is.null(fit))NA_integer_ else fit$converged, message = paste(unique(notes), collapse = '; '))
      }
      checks <- do.call(rbind, checks)
      good <- which(is.finite(checks$stress))
      checks$selected <- FALSE
      if (length(good)) {
        best <- good[which.min(checks$stress[good])]
        checks$selected[best] <- TRUE
        pts <- fits[[best]]$points[rownames(x),, drop = FALSE]
        scores[[id]] <- data.frame(stage, method, record_id = rownames(x), forest = group, NMDS1 = pts[, 1], NMDS2 = pts[, 2], stress = checks$stress[best])
      }
      diagnostics[[id]] <- checks
    }
  }
  distance_summary <- do.call(rbind, summaries)
  availability <- do.call(rbind, availability)
  taxa <- do.call(rbind, taxa);
  scores <- do.call(rbind, scores)
  family_summary <- do.call(rbind, families)
  table_s1b <- family_summary[order(match(family_summary$stage, stages),
      match(family_summary$forest, forests), family_summary$family),]
  put(data.frame(
      'Size class' = table_s1b$stage,
      'Forest' = table_s1b$forest,
      'Family' = table_s1b$family,
      'Taxa' = table_s1b$recorded_taxa,
      'Assigned taxa' = table_s1b$assigned_taxa_denominator,
      'Taxon share (%)' = ifelse(is.na(table_s1b$taxon_share_percent), NA_character_,
        sprintf('%.2f', table_s1b$taxon_share_percent)),
      'Stems' = table_s1b$recorded_stems,
      'Assigned stems' = table_s1b$assigned_stems_denominator,
      'Stem share (%)' = ifelse(is.na(table_s1b$stem_share_percent), NA_character_,
        sprintf('%.2f', table_s1b$stem_share_percent)),
      check.names = FALSE
    ), 'Table_S2_Family_representation_by_fragment_and_size_class')
  assignment_audit <- do.call(rbind, family_audit)
  put(assignment_audit, 'family_assignment_audit')
  write.csv(assignment_audit, file.path(derived, 'family_assignment_audit.csv'), row.names = FALSE, na = 'NA')
  put(do.call(rbind, richness), 'observed_taxon_richness')
  shared_taxa <- do.call(rbind, occurrences)
  shared_taxa$taxon <- sub('^Ocotea usambarensis$', 'Kuloa usambarensis',
    shared_taxa$taxon, ignore.case = TRUE)
  put(shared_taxa, 'taxa_shared_and_forest_only')
  put(availability, 'stage_record_summary');
  put(do.call(rbind, pairs), 'pairwise_dissimilarities')
  table_1 <- expand.grid(
    method = c('Jaccard', 'Bray-Curtis'),
    stage = c('Recorded juveniles', 'Listed larger stems'),
    stringsAsFactors = FALSE
  )
  comparisons <- c('Within Ngangao', 'Within Mbololo', 'Between fragments')
  formatted <- lapply(comparisons, function(comparison) {
      d <- distance_summary[which(distance_summary$comparison == comparison),]
      d <- d[match(paste(table_1$stage, table_1$method), paste(d$stage, d$method)),]
      ifelse(is.finite(d$median), sprintf('%.2f [%.2f, %.2f]', d$median, d$q1, d$q3), NA_character_)
    })
  put(data.frame(
      'Size class' = rep(c('Juveniles', 'Larger stems'), each = 2L),
      'Measure' = rep(c('Jaccard', 'Bray-Curtis'), 2L),
      'Within Ngangao' = formatted[[1L]],
      'Within Mbololo' = formatted[[2L]],
      'Between forests' = formatted[[3L]],
      check.names = FALSE
    ), 'Table_S1_Pairwise_dissimilarity_by_size_class', fileEncoding = 'UTF-8')
  table_s1 <- taxa
  table_s1$taxon <- paste0(toupper(substr(taxa$taxon, 1, 1)), tolower(substring(taxa$taxon, 2)))
  table_s1$taxon <- sub('^Ocotea usambarensis$', 'Kuloa usambarensis', table_s1$taxon)
  value_columns <- c('Ngangao_percent', 'Mbololo_percent', 'difference_pp')
  table_s1[value_columns] <- lapply(table_s1[value_columns], function(values)
    ifelse(is.na(values), NA_character_, sprintf('%.2f', values)))
  put(table_s1, 'Table_S3_Mean_within_record_relative_abundance_of_woody_taxa')
  put(do.call(rbind, diagnostics), 'nmds_diagnostics')
  if (!is.null(scores))put(scores, 'nmds_scores')
  save_figure <- function(name, draw, w, h, pointsize = 12)for (ext in c('pdf', 'png')) {
    path <- rq_file(out, paste0(name, '.', ext))
    if (ext == 'pdf')pdf(path, width = w, height = h, pointsize = pointsize, useDingbats = FALSE)else png(path, width = w * 200, height = h * 200, res = 200, pointsize = pointsize)
    tryCatch(draw(), finally = dev.off())
  }
  ordination <- function(method) {
    par(mfrow = c(1, 2), mar = c(4.4, 4.4, 5.4, 1), oma = c(2, 0, 2, 0))
    for (k in seq_along(stages)) {
      z <- if (!is.null(scores))scores[scores$stage == stages[k] & scores$method == method,]else NULL
      title <- paste0(LETTERS[k], '. ', stages[k])
      if (is.null(z) || !nrow(z)) {plot.new();
        title(main = title);
        text(.5, .5, 'NMDS unavailable; see diagnostics');
        next}
      plot(z$NMDS1, z$NMDS2, asp = 1, xlab = 'NMDS1', ylab = 'NMDS2', pch = symbols[z$forest],
        col = cols[z$forest], main = title, cex = 1.3)
      mtext(sprintf('Stress = %.3f', z$stress[1]), line = .35, cex = .9)
      mtext(paste(paste0(forests, ' n = ', as.integer(table(factor(z$forest, levels = forests)))), collapse = '; '),
        line = 1.5, cex = .8)
    }
    mtext(paste(method, 'NMDS of surveyed neighbourhoods'), side = 3, outer = TRUE, font = 2)
    mtext('Tan circles: Ngangao. Teal triangles: Mbololo. Axes are fitted separately by size class.', side = 1, outer = TRUE, cex = .85)
  }
  save_figure('Figure_1_Bray_Curtis_NMDS', function()ordination('Bray-Curtis'), 11, 5.5)
  save_figure('Figure_S1_Jaccard_NMDS', function()ordination('Jaccard'), 11, 5.5)
  save_figure('Figure_S5_Taxa_with_the_largest_relative_abundance_differences', function() {
      fills <- cols
      fills['Mbololo'] <- adjustcolor(cols['Mbololo'], alpha.f = .65)
      par(mfrow = c(1, 2), mar = c(3.9, 11.8, 2.3, .55), oma = c(0, 0, 0, 0),
        mgp = c(2.35, .55, 0), cex = 1, cex.axis = 1, cex.lab = 1)
      shown <- taxa[taxa$displayed,]
      xmax <- max(shown$Ngangao_percent, shown$Mbololo_percent) + 3
      for (k in seq_along(stages)) {
        par(mar = c(3.9, 11.8, 2.3, if (k == 1L) 1.8 else .55))
        if (k == 1L) par(mai = par('mai') + c(0, 4.5 / 25.4, 0, 0))
        d <- shown[shown$stage == stages[k],];
        y <- rev(seq_len(nrow(d)))
        x_limits <- if (k == 2L)c( - 1.5, 20)else c( - 3, xmax)
        plot(d$Ngangao_percent, y, type = 'n', xlim = x_limits, ylim = c(.65, nrow(d) + .9), yaxt = 'n',
          xaxt = if (k == 2L)'n'else 's', xaxs = 'i', yaxs = 'i',
          xlab = 'Mean within-record relative abundance (%)', ylab = '', cex.axis = 14 / 12, cex.lab = 18 / 12, font.lab = 2)
        bounds <- par('usr')
        strip_top <- grconvertY(grconvertY(bounds[4], from = 'user', to = 'inches') + .29,
          from = 'inches', to = 'user')
        rect(bounds[1], bounds[4], bounds[2], strip_top, col = 'grey90', border = 'black', xpd = NA)
        text(bounds[1], mean(c(bounds[4], strip_top)), paste0('(', letters[k], ') ', stages[k]),
          pos = 4, offset = .35, font = 2, cex = 17 / 12, xpd = NA)
        if (k == 2L)axis(1, at = seq(0, 20, 5), cex.axis = 14 / 12)
        labels <- as.expression(lapply(d$taxon, function(taxon) {
              name <- paste0(toupper(substr(taxon, 1, 1)), substring(taxon, 2))
              name <- sub('^Ocotea usambarensis$', 'Kuloa usambarensis', name, ignore.case = TRUE)
              bquote(italic(.(name)))
            }))
        axis(2, at = y, labels = labels, las = 2, cex.axis = 15.7 / 12)
        segments(d$Ngangao_percent, y, d$Mbololo_percent, y, col = 'grey70')
        points(d$Ngangao_percent, y, pch = 21, bg = cols[1], col = 'black', cex = 5.3, xpd = k == 2L)
        points(d$Mbololo_percent, y, pch = 24, bg = fills[2], col = 'black', cex = 5.6, xpd = k == 2L)
        text(d$Ngangao_percent, y, sprintf('%.1f', d$Ngangao_percent), cex = 14 / 12, font = 2, xpd = k == 2L)
        text(d$Mbololo_percent, y - .07, sprintf('%.1f', d$Mbololo_percent), cex = 14 / 12, font = 2, xpd = k == 2L)
        differences <- sprintf('[%+.1f]', d$difference_pp)
        difference_cex <- 14 / 12
        ngangao_x <- grconvertX(d$Ngangao_percent, from = 'user', to = 'inches')
        mbololo_x <- grconvertX(d$Mbololo_percent, from = 'user', to = 'inches')
        # Measure clearance from the marker outlines in device inches.
        midpoint <- (ngangao_x + mbololo_x) / 2
        fits <- abs(ngangao_x - mbololo_x) >= strwidth(differences, units = 'inches', cex = difference_cex) + 2 * .28
        right_edge <- pmax(ngangao_x + .20, mbololo_x + .28)
        label_x <- grconvertX(ifelse(fits, midpoint, right_edge + 2 / 25.4), from = 'inches', to = 'user')
        label_y <- grconvertY(grconvertY(y, from = 'user', to = 'inches') + .055, from = 'inches', to = 'user')
        if (any(fits))
          text(label_x[fits], label_y[fits], differences[fits], adj = c(.5, 0), cex = difference_cex)
        if (any(!fits))
          text(label_x[!fits], y[!fits], differences[!fits], adj = c(0, .5), cex = difference_cex, xpd = TRUE)
        legend(bounds[2] - .02 * diff(bounds[1:2]), 2, forests, xjust = 1, yjust = .5,
          pch = c(21, 24), pt.bg = fills, col = 'black', pt.cex = 3,
          bty = 'o', box.col = 'black', cex = 15.5 / 12, y.intersp = 1.75, x.intersp = 1.2)
      }
    }, 15.4, 9.12)
  asset_dir <- file.path(root, 'analysis/Datasets/figure_assets/phylopic')
  credits <- read.csv(file.path(asset_dir, 'credits.csv'), stringsAsFactors = FALSE, fileEncoding = 'UTF-8')
  stopifnot(!anyDuplicated(credits$family))
  icons <- setNames(lapply(file.path(asset_dir, credits$file), png::readPNG, native = TRUE), credits$family)
  leading_families <- lapply(stages, function(stage) {
      d <- family_summary[family_summary$stage == stage & is.finite(family_summary$taxon_share_percent),]
      rev(head(names(sort(tapply(d$recorded_taxa, d$family, sum), decreasing = TRUE)), 8))
    })
  shown_families <- sort(unique(unlist(leading_families)))
  image_credits <- merge(data.frame(family = shown_families), credits, by = 'family', all.x = TRUE, sort = TRUE)
  image_credits$status <- ifelse(is.na(image_credits$file), 'No matching PhyloPic silhouette located',
    'Representative of the family')
  put(image_credits, 'Figure_1_PhyloPic_credits')
  save_figure('Figure_2_NMDS_of_juvenile_and_listed_larger_stem_assemblages', function() {
      par(mfrow = c(2, 2), mar = c(3.5, 4.1, 2.2, .65), oma = c(0, 0, 0, 0),
        mgp = c(2.45, .8, 0), cex = 1, cex.axis = 15 / 12, cex.lab = 15 / 12, font.lab = 2)
      methods <- c('Jaccard', 'Bray-Curtis')
      for (m in seq_along(methods)) {
        method <- methods[m]
        for (k in seq_along(stages)) {
          stage <- stages[k]
          z <- if (!is.null(scores))scores[scores$stage == stage & scores$method == method,]else NULL
          heading <- paste0('(', letters[k + 2L * (m - 1L)], ') ', stage)
          if (!is.null(z) && nrow(z)) {
            y_limits <- range(z$NMDS2)
            y_limits <- y_limits + c( - .22, .15) * diff(y_limits)
            x_limits <- range(z$NMDS1)
            x_limits <- x_limits + c( - .04, .025) * diff(x_limits)
            plot(z$NMDS1, z$NMDS2, asp = 1, col = 'black', bg = cols[z$forest],
              pch = c(Ngangao = 21, Mbololo = 24)[z$forest],
              xlim = x_limits, ylim = y_limits, xaxs = 'i', xlab = 'NMDS1', ylab = 'NMDS2', cex = 2, lwd = 1)
            title(main = heading, adj = 0, cex.main = 16 / 12)
            stress_label <- bquote(.(paste0(method, ' NMDS stress = ')) *
              bold(.(sprintf('%.3f', z$stress[1]))))
            legend('topleft', as.expression(stress_label),
              bty = 'n', inset = .025, cex = 15 / 12)
            legend('bottomleft', forests, col = 'black', pt.bg = cols[forests],
              pch = c(Ngangao = 21, Mbololo = 24)[forests], bty = 'n', pt.cex = 1.6, inset = .025, cex = 15 / 12)
          } else {
            plot.new();
            title(main = heading, adj = 0, cex.main = 16 / 12);
            text(.5, .5, 'NMDS unavailable', cex = 15 / 12)
          }
        }
      }
    }, 15.18, 8.92584)
  save_figure('Figure_S4_Family_representation_in_Ngangao_and_Mbololo', function() {
      par(mfrow = c(1, 2), omi = c(.45, 0, 0, 0), mgp = c(2.45, .8, 0),
        cex = 1, cex.axis = 15.5 / 12, font.lab = 2)
      for (k in seq_along(stages)) {
        stage <- stages[k]
        heading <- paste0('(', letters[k], ') ', stage)
        d <- family_summary[family_summary$stage == stage & is.finite(family_summary$taxon_share_percent),]
        shown <- leading_families[[k]]
        par(mar = c(3.6, 10.3, 2.2, .85))
        if (!length(shown)) {plot.new();
          title(main = heading, adj = 0, cex.main = 16.5 / 12);
          text(.5, .5, 'No assigned families', cex = 15.5 / 12);
          next}
        y <- seq_along(shown)
        # Measure the silhouette-to-label gap on the output device in millimetres.
        plot(NA, xlim = c(0, 21), ylim = c(.48, max(y) + .52),
          axes = FALSE, xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
        title(main = heading, adj = 0, cex.main = 16.5 / 12)
        ticks <- seq(0, 21, 3)
        axis(1, at = ticks)
        label_right <- grconvertX(grconvertX(0, from = 'user', to = 'inches') - .06,
          from = 'inches', to = 'user')
        text(label_right, y, shown, adj = c(1, .5), xpd = NA, cex = 15.5 / 12)
        abline(v = ticks, col = '#EBEBEB', lwd = .7)
        abline(h = y, col = '#EBEBEB', lwd = .7)
        box(col = '#333333')
        for (j in seq_along(forests)) {
          f <- d[d$forest == forests[j],];
          f <- f[match(shown, f$family),]
          yy <- y + c(-.20, .20)[j]
          positive <- which(f$taxon_share_percent>0)
          rect(0, yy[positive] - .18, f$taxon_share_percent[positive], yy[positive] + .18,
            col = cols[j], border = 'black')
          labels <- sprintf('%.1f%%', f$taxon_share_percent)
          inside <- f$taxon_share_percent + .28 + strwidth(labels, font = 2, cex = 15.5 / 12)>21
          text(f$taxon_share_percent + ifelse(inside, - .28, .28), yy, labels,
            col = 'black',
            pos = ifelse(inside, 2, 4), offset = 0, font = 2, cex = 15.5 / 12)
        }
        # Give each representative an equally sized box and preserve its aspect ratio.
        units_per_inch <- c(diff(par('usr')[1:2]), diff(par('usr')[3:4])) / par('pin')
        for (j in seq_along(shown)) {
          icon <- icons[[shown[j]]]
          label_width <- strwidth(shown[j], units = 'inches', cex = 15.5 / 12)
          icon_right <- label_right - (label_width + 3 / 25.4) * units_per_inch[1]
          if (is.null(icon)) {
            text(icon_right - .15 * units_per_inch[1], y[j], '-', col = 'grey55', xpd = NA, cex = 15.5 / 12)
          } else {
            dimensions <- rev(dim(icon))
            radius <- (.18 + 1.5 / 144) * dimensions / max(dimensions) * units_per_inch
            rasterImage(icon, icon_right - 2 * radius[1], y[j] - radius[2],
              icon_right, y[j] + radius[2], xpd = NA)
          }
        }
      }
      text(grconvertX(.5, from = 'ndc', to = 'user'),
        grconvertY(.078, from = 'ndc', to = 'user'),
        'Share of family-assigned taxa (%)', font = 2, cex = 16.5 / 12, xpd = NA)
      legend(grconvertX(.5, from = 'ndc', to = 'user'),
        grconvertY(.032, from = 'ndc', to = 'user'),
        rev(forests), fill = cols[rev(forests)], border = 'black',
        horiz = TRUE, xjust = .5, yjust = .5, bty = 'n', cex = 17.5 / 12, xpd = NA)
    }, 15, 7)
  invisible(list(output = out, comparison = distance_summary, provisional = provisional))
}

# ============================================================
# 8. Run floristic composition
# ============================================================

run_floristic_composition <- function(root, here = file.path(root, 'analysis/scripts'), repeats = 8L, refresh = FALSE) {
  floristic_temp_ready()
  local_lib <- file.path(here, '.r-library')
  if (dir.exists(local_lib)).libPaths(c(local_lib, .libPaths()))
  for (pkg in c('readxl', 'dplyr', 'vegan', 'rWCVP', 'rWCVPdata', 'png'))if (!requireNamespace(pkg, quietly = TRUE))stop('Missing R package: ', pkg)
  out <- rq_dirs(root, 1)
  inputs <- floristic_inputs(root, out)
  raw <- floristic_provisional(root, out, export = FALSE)
  if (is.null(inputs)) {
    write.csv(raw$records, rq_file(out, 'worksheet_selection.csv'), row.names = FALSE)
    write.csv(raw$source_rows, rq_file(out, 'source_rows.csv'), row.names = FALSE)
    write.csv(raw$evidence, rq_file(out, 'worksheet_stage_evidence_audit.csv'), row.names = FALSE)
    return(floristic_composition_report(raw$records, raw$source_rows, raw$evidence, out, TRUE, repeats, root = root))
  }
  rec <- data.frame(record_id = inputs$q$plot_id, forest = inputs$q$forest, included = inputs$q$primary_eligible)
  floristic_composition_report(rec, raw$source_rows, raw$evidence, out, FALSE, repeats, list(inputs$J, inputs$A), root = root)
}

# ============================================================
# 9. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep('^--file=', commandArgs(), value = TRUE)
  start <- if (length(arg))dirname(normalizePath(sub('^--file=', '', arg[1])))else getwd()
  root <- floristic_project_root(start)
  here <- file.path(root, 'analysis/scripts')
  run_floristic_composition(root, here)
  source(file.path(file.path(root, 'analysis'), 'scripts', 'combine_output_tables.R'), local = TRUE)
  combine_output_tables(file.path(file.path(root, 'analysis'), 'outputs', 'research question 1'), '01_woody_composition.R')
}
