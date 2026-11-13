# ============================================================
# 0. Sampling design
# ============================================================

# Run with Rscript analysis/scripts/sampling_design.R.
sampling_design <- function(root) {
  out <- file.path(root, "analysis", "outputs", "sampling_design", "figures")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  metres <- seq(0, 400, by = 50)
  side <- c(1, - 1, - 1, 1, - 1, 1, 1, - 1, 1)
  assets <- file.path(root, "analysis", "Datasets", "figure_assets", "phylopic")
  images <- file.path(assets, c("sampling_conifer.png", "Rubiaceae.png"))
  source_urls <- paste0("https://www.phylopic.org/images/",
    c("1fe9c2ad-6083-4829-acd0-e4a9e67608d1", "c2c6dc22-04e1-41f2-8137-2b995ef1d8ab"))
  representatives <- c("Podocarpus", "Coffea alleizettei")
  stopifnot(requireNamespace("png", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
  artwork <- lapply(images, function(path) {
      img <- png::readPNG(path)
      rgb <- grDevices::col2rgb("#245D4A") / 255
      alpha <- img[, , dim(img)[3]]
      array(c(rep(rgb, each = length(alpha)), alpha), c(dim(alpha), 4L))
    })
  draw <- function() {
    layout(matrix(1:3, ncol = 1), heights = c(1.15, 1.8, 1.15))
    par(mar = c(.6, 1.1, 1.8, .2), family = "sans", xaxs = "i", yaxs = "i")
    blue <- "#175777"
    green <- "#23654F"
    grey <- "#4A5561"
    font_increase <- 2 / 13
    panel <- function(title, xlim, ylim, cex = 1.45) {
      plot.new()
      plot.window(xlim, ylim)
      title(main = title, adj = 0, line = .35, cex.main = cex, col.main = "black")
    }
    box <- function(x, y, w, h, label, fill, border = fill, cex = 1.22, text_colour = "#202B34") {
      rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = border)
      text(x, y, label, cex = cex, col = text_colour)
    }
    panel("(a) Focal species searches", c(0, 100), c(0, 25))
    box(22, 18.5, 40, 11, "Initial searches in\nNgangao and Mbololo", "#E6F0F6", cex = 1.4 + font_increase)
    box(79, 18.5, 40, 11, "Targeted visits guided by\nlocal community members", "#E6F0F6", cex = 1.4 + font_increase)
    arrows(43, 18.5, 58, 18.5, length = .09, col = blue)
    badge <- function(x, i) {
      ratio <- (diff(par("usr")[1:2]) / par("pin")[1]) /
      (diff(par("usr")[3:4]) / par("pin")[2])
      h <- 12.4
      w <- h * ncol(artwork[[i]]) / nrow(artwork[[i]]) * ratio
      symbols(x, 6.5, circles = .30, inches = .30, bg = "#E6F1E9", fg = NA, add = TRUE)
      rasterImage(artwork[[i]], x - w / 2, 6.5 - h / 2, x + w / 2, 6.5 + h / 2)
    }
    badge(24, 1)
    badge(62, 2)
    text(28, 6.5, expression(italic(Afrocarpus ~ usambarensis)), adj = 0, cex = 1.25)
    text(66, 6.5, expression(italic(Coffea ~ fadenii)), adj = 0, cex = 1.25)
    panel("(b) Intended transect and quadrat layout", c( - 30, 455), c( - 46, 56), cex = 1.45 + font_increase)
    rect( - 20, - 25, 420, 25, col = "#F2F7F3", border = NA)
    segments(0, 0, 400, 0, col = blue, lwd = 2.5)
    segments(metres, 0, metres, 12 * side, col = grey, lty = 2)
    rect(metres - 5, pmin(12 * side, 22 * side), metres + 5,
      pmax(12 * side, 22 * side), col = "#B9DAC9", border = green)
    points(metres, rep(0, length(metres)), pch = 19, col = blue)
    text(metres, rep( - 29, length(metres)), metres, cex = 1.12 + font_increase)
    text(0, - 40, "Start", cex = 1.5 + font_increase, font = 2, col = "black")
    text(426, - 29, "m", cex = 1.12 + font_increase)
    arrows(0, 40, 400, 40, code = 3, angle = 25, length = .07, col = blue)
    text(200, 48, "400 m transect", font = 2, cex = 1.25 + font_increase)
    text(200, - 40, "Nine points at 50 m intervals, with one 10 m by 10 m quadrat per point", cex = 1.20 + font_increase)
    text(0, 30, "Quadrat sides shown as an example of random assignment", adj = 0, cex = 1.2 + font_increase, col = "black")
    panel("(c) Survey plots used in the comparison", c(0, 100), c(18, 55), cex = 1.45 + font_increase)
    box(50, 48, 30, 11, "28 plots", "#E6F0F6", cex = 1.4 + font_increase, text_colour = "black")
    segments(50, 42.5, 50, 40, col = grey)
    segments(25, 40, 75, 40, col = grey)
    arrows(c(25, 75), 40, c(25, 75), 36, length = .07, col = grey)
    box(25, 30, 35, 12, "Ngangao\n9 plots, all retained", "#D4B58A", cex = 1.4 + font_increase, text_colour = "black")
    box(75, 28, 44, 19, "Mbololo\n19 plots, 14 retained and 5 excluded\n(roadside or farm plots)", "#508A92", cex = 1.4 + font_increase, text_colour = "black")
  }
  stem <- file.path(out, "Figure_S2_Field_sampling_in_Ngangao_and_Mbololo")
  grDevices::png(paste0(stem, ".png"), width = 3300, height = 1872, res = 300, pointsize = 13)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 11, height = 6.24, pointsize = 13)
  tryCatch(draw(), finally = grDevices::dev.off())
  silhouette <- function(i) {
    img <- readBin(images[i], "raw", n = file.info(images[i])$size)
    sprintf('<img src="data:image/png;base64,%s" alt="PhyloPic silhouette of %s">',
      jsonlite::base64_enc(img), representatives[i])
  }
  stations <- vapply(seq_along(metres), function(i) {
      sprintf(paste0('<button class="station %s" type="button" aria-pressed="%s" ',
          'data-station="%d" data-distance="%d" data-side="%s"><span class="quadrat"></span>',
          '<span class="point"></span><span class="distance">%d m</span></button>'),
        if (side[i] > 0) "upper" else "lower", if (i == 1L) "true" else "false",
        i, metres[i], if (side[i] > 0) "above" else "below", metres[i])
    }, character(1))
  writeLines(c(
      '<!doctype html><html lang="en"><meta charset="utf-8">',
      '<meta name="viewport" content="width=device-width,initial-scale=1">',
      '<title>Taita Hills sampling design</title>',
    '<style>
    *{box-sizing:border-box}body{margin:0;background:#edf3ef;color:#19362d;font:18px/1.5 system-ui,sans-serif}
    main{max-width:960px;margin:auto;padding:clamp(14px,4vw,36px);padding-bottom:12px}
    header{background:#173c31;color:#fff;padding:24px 30px;border-radius:18px}
    h1,h2,p{margin:0 0 10px}h1{font-size:clamp(25px,5vw,40px)}h2{font-size:24px}
    section{background:#fff;border:1px solid #d3e1d7;border-radius:16px;margin-top:16px;padding:20px}
    .cards{display:flex;flex-wrap:wrap;gap:12px}.card{flex:1;min-width:190px;padding:12px;border-radius:10px;background:#eaf3ec;display:flex;gap:12px;align-items:center}
    .card img{width:110px;height:110px;object-fit:contain;flex:none;padding:8px;background:#d9eadf;border-radius:50%;filter:sepia(.32)}
    .track{display:grid;grid-template-columns:repeat(9,1fr);position:relative;min-height:134px}
    .track:before{content:"";position:absolute;top:56px;left:5.55%;right:5.55%;border-top:2px solid #28654d}
    .station{position:relative;min-width:0;border:0;background:transparent;cursor:pointer;color:#173c31;font:inherit;padding:0}
    .station:focus-visible{outline:3px solid #175777;outline-offset:-2px;border-radius:6px}
    .point,.quadrat,.distance{position:absolute;left:50%;transform:translateX(-50%)}
    .point{top:50px;border-radius:50%;width:14px;height:14px;background:#28654d;border:2px solid white}
    .quadrat{width:21px;height:21px;border:1px solid #28654d;background:#b9dac9;top:16px}
    .lower .quadrat{top:76px}.distance{bottom:8px;font-size:clamp(10px,2.5vw,17px);white-space:nowrap}
    .station:hover .quadrat,.station[aria-pressed=true] .quadrat{background:#f1cb78;border:2px solid #875713}
    .station[aria-pressed=true] .point{background:#875713}
    output{display:block;background:#f6f1e1;padding:12px;border-radius:8px}
    details{border-top:1px solid #d3e1d7;padding:10px 0}summary{cursor:pointer;font-weight:650}
    .note{font-size:16px;color:#41564d}
    </style>',
      '<main><header><h1>Taita Hills sampling design</h1><p>Ngangao and Mbololo forest fragments</p></header>',
      '<section><h2>(a) Focal species searches</h2><p>Local community members guided targeted visits after the initial searches for focal plants and associated woody vegetation.</p><div class="cards">',
      paste0('<div class="card">', silhouette(1), '<i>Afrocarpus usambarensis</i></div>'),
      paste0('<div class="card">', silhouette(2), '<i>Coffea fadenii</i></div>'),
      '<section><h2>(b) Intended transect and quadrat layout</h2><p>Nine stations at 50 m intervals along a 400 m transect. Select a station to inspect its example quadrat.</p>',
      '<div class="track" role="group" aria-label="Transect sampling stations">', stations, '</div>',
      '<output id="station-info" aria-live="polite">Station 1 at 0 m. Example 10 m by 10 m quadrat above the transect line.</output>',
      '<p class="note">Quadrat sides are illustrative. The locations along the transect, surveyed footprints, overlap and independence could not be verified.</p></section>',
      '<section><h2>(c) Worksheet screening</h2><p>23 of 28 plots entered the primary descriptive comparison. Counts were reported per worksheet.</p>',
      '<details><summary>Ngangao &middot; 9 of 9 retained</summary>All nine plots entered the primary comparison.</details>',
      '<details><summary>Mbololo &middot; 14 of 19 retained</summary>Two roadside or farm plots and three plots with unresolved identities were excluded.</details>',
      '</section></main>',
    '<script>
    const buttons = document.querySelectorAll(".station");
    const info = document.getElementById("station-info");
    buttons.forEach(button => button.addEventListener("click", () => {
      buttons.forEach(item => item.setAttribute("aria-pressed", String(item === button)));
      const data = button.dataset;
      info.textContent = "Station " + data.station + " at " + data.distance +
        " m. Example 10 m by 10 m quadrat " + data.side + " the transect line.";
    }));
    </script></html>'
    ), paste0(stem, ".html"), useBytes = TRUE)
  invisible(stem)
}

# ============================================================
# 1. Run the script
# ============================================================

if (sys.nframe() == 0L) {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  start <- if (length(arg)) dirname(sub("^--file=", "", arg[1])) else getwd()
  root <- normalizePath(start, winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(root, "analysis/scripts/sampling_design.R"))) {
    parent <- dirname(root)
    if (identical(root, parent)) stop("Cannot locate the Taita project.")
    root <- parent
  }
  sampling_design(root)
}
