# ============================================================
# 0. Data collection
# ============================================================

# Run from the project with Rscript analysis/scripts/data_collection.R.
# Read the original field workbooks and write one auditable data workbook.
data_collection <- function(start = getwd()) {
  root <- normalizePath(start, winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(root, "analysis/scripts/01_woody_composition.R"))) {
    parent <- dirname(root)
    if (identical(parent, root)) stop("Cannot locate the Taita project.")
    root <- parent
  }
  parser <- new.env()
  sys.source(file.path(root, "analysis/scripts/01_woody_composition.R"), envir = parser)
  input <- parser$floristic_provisional(root, out = NULL, export = FALSE)
  plots <- input$records
  rows <- input$source_rows
  stems <- rows[rows$stage == "Listed larger stems", ]
  juveniles <- rows[rows$stage == "Recorded juveniles", ]
  focal_taxa <- c("coffea fadenii", "afrocarpus usambarensis")
  number <- function(x, unit = "") {
    x <- sub(paste0("[[:space:]]*", unit, "$"), "", trimws(x), ignore.case = TRUE)
    suppressWarnings(as.numeric(sub(",", ".", x, fixed = TRUE)))
  }
  degrees <- function(x, sign) {
    tokens <- regmatches(x, gregexpr("[0-9]+(?:[.,][0-9]+)?", x, perl = TRUE))
    vapply(tokens, function(token) {
        if (length(token) != 2L) return(NA_real_)
        value <- number(token)
        if (value[2] >= 60) return(NA_real_)
        sign * (value[1] + value[2] / 60)
      }, numeric(1))
  }
  stems$dbh_cm <- number(stems$raw_dbh, "cm")
  stems$height_m <- number(stems$raw_height, "m")
  stems$dbh_above_5_cm <- ifelse(is.finite(stems$dbh_cm), stems$dbh_cm > 5, NA)
  stems$is_focal <- stems$recorded_label %in% focal_taxa
  juveniles$is_focal <- juveniles$recorded_label %in% focal_taxa
  stems$inventory_number <- rep(NA_character_, nrow(stems))
  juveniles$inventory_number <- juveniles$raw_count <- rep(NA_character_, nrow(juveniles))
  plots$latitude_dd <- degrees(plots$latitude_label, - 1)
  plots$longitude_dd <- degrees(plots$longitude_label, 1)
  cover_labels <- c(canopy = "^Upper canopy cover", shrub = "^Cover of shrub layer",
    herbaceous = "^Cover of herbaceous layer", litter = "^Cover of litter",
    bare_soil_rocks = "^Cover of bare")
  gps <- list()
  headers <- list()
  cover <- list()
  for (i in seq_len(nrow(plots))) {
    plot <- plots[i, ]
    sheet <- readxl::read_excel(file.path(root, plot$workbook), sheet = plot$sheet,
      range = readxl::cell_limits(c(1, 1), c(NA, 16)), col_names = FALSE,
      col_types = "text", .name_repair = "minimal")
    header_cells <- trimws(unlist(sheet[seq_len(min(10L, nrow(sheet))), ], use.names = FALSE))
    field <- function(pattern) {
      values <- header_cells[which(grepl(pattern, header_cells, ignore.case = TRUE))]
      if (length(values)) paste(values, collapse = " | ") else NA_character_
    }
    for (layer in names(cover_labels)) {
      raw <- field(cover_labels[[layer]])
      value <- number(sub("^[^:]*:[[:space:]]*", "", raw), "%")
      cover[[paste(i, layer)]] <- data.frame(record_id = plot$record_id, layer,
        cover_percent = value, raw_field = raw,
        status = if (is.na(raw)) "Blank" else if (!is.finite(value)) "Review raw field"
        else if (value < 0 || value > 100) "Outside 0–100%" else "Recorded")
    }
    stem_idx <- which(stems$record_id == plot$record_id)
    juvenile_idx <- which(juveniles$record_id == plot$record_id)
    stems$inventory_number[stem_idx] <- sheet[[1]][stems$source_row[stem_idx]]
    juveniles$inventory_number[juvenile_idx] <- sheet[[6]][juveniles$source_row[juvenile_idx]]
    juveniles$raw_count[juvenile_idx] <- sheet[[8]][juveniles$source_row[juvenile_idx]]
    taxon <- plot$focal_header_taxon
    raw_taxon <- sub("^Name of the species:[[:space:]]*", "",
      field("^Name of the species:"), ignore.case = TRUE)
    raw_height <- sub("^Height:[[:space:]]*", "", field("^Height:"), ignore.case = TRUE)
    if (any(!is.na(c(raw_taxon, plot$focal_header_dbh, raw_height,
            plot$focal_header_saplings, plot$focal_header_larger_count)))) {
      dbh <- number(plot$focal_header_dbh, "cm")
      listed <- stems[stem_idx, ]
      matches <- listed$source_row[which(listed$recorded_label == taxon &
          is.finite(listed$dbh_cm) & abs(listed$dbh_cm - dbh) < 1e-8)]
      inventory_header <- input$evidence$header_row[
        input$evidence$record_id == plot$record_id & input$evidence$stage == "Listed larger stems"]
      status <- if (is.na(taxon) || !nzchar(taxon)) "Header taxon unavailable"
      else if (!is.finite(dbh) || dbh <= 0) "No usable header DBH"
      else if (is.na(inventory_header)) "Numbered inventory unavailable"
      else if (length(matches)) "Possible duplicate"
      else "No same-taxon DBH match"
      headers[[i]] <- data.frame(record_id = plot$record_id, raw_taxon, normalized_taxon = taxon,
        raw_dbh = plot$focal_header_dbh, dbh_cm = dbh, raw_height,
        height_m = number(raw_height, "m"), raw_sapling_count = plot$focal_header_saplings,
        raw_mature_count = plot$focal_header_larger_count, duplicate_status = status,
        matching_worksheet_rows = paste(matches, collapse = ", "))
    }
    gps_header <- which(tolower(trimws(sheet[[10]])) == "no" &
      tolower(trimws(sheet[[11]])) %in% c("lat", "latitude") &
      tolower(trimws(sheet[[12]])) %in% c("long", "longitude"))
    if (length(gps_header) > 1L) stop("Ambiguous GPS table in ", plot$record_id)
    if (length(gps_header)) {
      gps_rows <- seq_len(nrow(sheet))
      gps_rows <- gps_rows[gps_rows > gps_header]
      end <- which(!grepl("^[0-9]+$", trimws(sheet[[10]][gps_rows])))
      if (length(end)) gps_rows <- head(gps_rows, end[1] - 1L)
      if (length(gps_rows)) {
        gps[[i]] <- data.frame(record_id = plot$record_id, source_row = gps_rows,
          observation_number = sheet[[10]][gps_rows], context = sheet[[10]][gps_header - 1L],
          header_taxon = taxon, raw_latitude = sheet[[11]][gps_rows],
          raw_longitude = sheet[[12]][gps_rows],
          latitude_dd = degrees(sheet[[11]][gps_rows], - 1),
          longitude_dd = degrees(sheet[[12]][gps_rows], 1),
          raw_elevation = sheet[[13]][gps_rows], raw_mature_count = sheet[[14]][gps_rows],
          raw_sapling_count = sheet[[15]][gps_rows], raw_dbh = sheet[[16]][gps_rows])
      }
    }
  }
  focal_counts <- do.call(rbind, lapply(seq_len(nrow(input$evidence)), function(i) {
        evidence <- input$evidence[i, ]
        stage_rows <- rows[rows$record_id == evidence$record_id & rows$stage == evidence$stage, ]
        do.call(rbind, lapply(focal_taxa, function(taxon) {
              selected <- stage_rows[stage_rows$recorded_label == taxon, ]
              count <- if (nrow(stage_rows) && all(stage_rows$usable_count)) sum(selected$count) else NA_real_
              data.frame(record_id = evidence$record_id, taxon, stage = evidence$stage,
                recorded_count = count, matching_rows = nrow(selected), table_status = evidence$status)
            }))
      }))
  notes <- data.frame(
    topic = c("Scope", "Plot coordinates", "Nearby GPS records", "Cover", "Numbered stems",
      "Recorded juveniles", "Juvenile heights", "Focal headers", "Duplicate check", "Blank cells"),
    interpretation = c(
      "All plot worksheets from the Ngangao and Mbololo field workbooks are represented. Plot IDs join the sheets.",
      "The Garmin eTrex 30 coordinates in the field headers locate corresponding plots. Decimal degrees use south and east; the datum was not independently verified.",
      "Separate numbered GPS observations around plots retain their worksheet context. Header taxon is provided for context, not assigned as an identification for each observation.",
      "Cover layers are separate percentage fields. Unparsed or out-of-range entries retain their original wording.",
      "Each numbered row represents a listed stem. The DBH >5 cm protocol is flagged; stems at or below 5 cm and missing measurements remain in the workbook.",
      "The numbered seedling/sapling tables contain combined counts by taxon. Focal counts use exact normalized field labels and exclude separate header and GPS entries.",
      "No dedicated juvenile height field was identified in the numbered seedling/sapling tables in these two workbooks. Header heights and numbered stem heights are separate.",
      "Header observations retain the recorded taxon, measurements and counts separately from the numbered inventory.",
      "A possible duplicate is an exact normalized taxon and matching numeric DBH on the same worksheet. A matching row does not establish individual identity.",
      "Blank cells indicate unavailable or unparsed values. A recorded zero is stored as zero."))
  sheets <- list("Plots" = plots[, c("record_id", "forest", "workbook", "sheet",
        "worksheet_plot_number", "header_plot_number", "header_forest", "survey_date_label",
        "latitude_label", "longitude_label", "latitude_dd", "longitude_dd", "worksheet_context")],
    "Focal counts" = focal_counts, "Cover" = do.call(rbind, cover),
    "Listed stems" = stems[, c("record_id", "source_row", "inventory_number", "source_label",
        "is_focal", "raw_dbh", "dbh_cm", "dbh_above_5_cm", "raw_height", "height_m")],
    "Recorded juveniles" = juveniles[, c("record_id", "source_row", "inventory_number",
        "source_label", "is_focal", "raw_count", "count", "usable_count")],
    "Focal headers" = do.call(rbind, headers), "Nearby GPS" = do.call(rbind, gps), "Notes" = notes)
  stopifnot(!anyDuplicated(plots$record_id), nrow(sheets$Cover) == 5L * nrow(plots),
    nrow(focal_counts) == nrow(plots) * 4L)
  out <- file.path(root, "analysis/outputs/data_collection")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out, "Taita_Hills_data_collection.xlsx")
  book <- openxlsx::createWorkbook()
  heading <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#245D4A",
    textDecoration = "bold", wrapText = TRUE, valign = "center")
  for (name in names(sheets)) {
    data <- sheets[[name]]
    openxlsx::addWorksheet(book, name, gridLines = FALSE)
    openxlsx::writeData(book, name, data, headerStyle = heading, withFilter = name != "Notes")
    openxlsx::freezePane(book, name, firstRow = TRUE, firstCol = name != "Notes")
    widths <- pmin(48, pmax(16, nchar(names(data)) + 4))
    openxlsx::setColWidths(book, name, cols = seq_along(widths), widths = widths)
    openxlsx::setRowHeights(book, name, rows = 1, heights = 30)
  }
  openxlsx::setColWidths(book, "Notes", cols = 2, widths = 95)
  openxlsx::addStyle(book, "Notes", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
    rows = 2:(nrow(notes) + 1L), cols = 2, gridExpand = TRUE)
  openxlsx::setRowHeights(book, "Notes", rows = 2:(nrow(notes) + 1L), heights = 34)
  openxlsx::saveWorkbook(book, path, overwrite = TRUE)
  message("Wrote ", path)
  invisible(path)
}

# ============================================================
# 1. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  data_collection(if (length(arg)) dirname(sub("^--file=", "", arg[1])) else getwd())
}
