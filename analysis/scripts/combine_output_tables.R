combine_output_tables <- function(out, script_name) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Install openxlsx to export the combined workbook.")
  }
  tables <- file.path(out, "tables")
  files <- list.files(tables, pattern = "[.]csv$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) return(invisible(NULL))
  files <- files[order(!grepl("^Table_", basename(files)), basename(files))]
  book <- openxlsx::createWorkbook()
  section <- openxlsx::createStyle(fgFill = "#214F56", fontColour = "#FFFFFF", textDecoration = "bold")
  header <- openxlsx::createStyle(fgFill = "#E3EFF0", textDecoration = "bold")
  for (i in seq_along(files)) {
    file <- files[i]
    label <- tools::file_path_sans_ext(basename(file))
    sheet <- paste0(sprintf("%02d ", i), substr(gsub("_", " ", label), 1L, 28L))
    openxlsx::addWorksheet(book, sheet, gridLines = FALSE)
    openxlsx::freezePane(book, sheet, firstActiveRow = 3)
    row <- 1L
    widths <- numeric()
    csv_text <- rawToChar(readBin(file, "raw", n = file.info(file)$size))
    if (validUTF8(csv_text)) {
      Encoding(csv_text) <- "UTF-8"
    } else {
      csv_text <- iconv(csv_text, from = "Windows-1252", to = "UTF-8")
    }
    data <- read.csv(text = csv_text, check.names = FALSE, colClasses = "character",
      na.strings = NULL)
    cols <- seq_along(data)
    if (!length(cols)) next
    openxlsx::writeData(book, sheet, basename(file), startRow = row, colNames = FALSE)
    openxlsx::addStyle(book, sheet, section, rows = row, cols = cols, gridExpand = TRUE)
    openxlsx::writeData(book, sheet, t(names(data)), startRow = row + 1L, colNames = FALSE)
    openxlsx::addStyle(book, sheet, header, rows = row + 1L, cols = cols, gridExpand = TRUE)
    if (nrow(data)) {
      openxlsx::writeData(book, sheet, unname(as.matrix(data)), startRow = row + 2L, colNames = FALSE)
    }
    for (j in cols) {
      widths[j] <- min(55, max(14, nchar(names(data)[j]) + 2,
        nchar(data[[j]], type = "width", allowNA = TRUE) + 2, na.rm = TRUE))
    }
    openxlsx::setColWidths(book, sheet, cols = seq_along(widths), widths = widths)
  }
  path <- file.path(tables, paste0(tools::file_path_sans_ext(script_name), "_combined_tables.xlsx"))
  if (!isTRUE(openxlsx::saveWorkbook(book, path, overwrite = TRUE, returnValue = TRUE))) {
    stop("Cannot save workbook. Close it in Excel and retry: ", path)
  }
  invisible(path)
}
