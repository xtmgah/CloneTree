#' @importFrom utils globalVariables
NULL

#' Default CloneTree palette
#'
#' A muted palette close to the published multi-region clone figures.
#'
#' @return A character vector of hex colors.
#' @export
clonetree_palette <- function() {
  c(
    "#8E8E8E", "#A72B30", "#80B88C", "#4F3668", "#CE834E",
    "#587560", "#68AFC3", "#E1E675", "#536D92", "#B44264",
    "#2A4070", "#5BA35E", "#A67353", "#DFBD52", "#7D839F",
    "#70482F", "#AB524F", "#C58AA8", "#6A9FB5", "#B7B46E"
  )
}

#' Use the CloneTree figure font
#'
#' Enables Roboto Condensed through `showtext` when the local font files are
#' available. If `showtext` or the font files are unavailable, CloneTree falls
#' back to the device default font.
#'
#' @param family Font family name.
#' @param regular,bold,italic,bolditalic Font file paths.
#' @param dpi DPI passed to `showtext::showtext_opts()`.
#'
#' @return Invisibly returns the active font family, or `""` for the device
#'   default.
#' @export
use_clonetree_fonts <- function(family = "Roboto Condensed",
                                regular = "/Users/zhangt8/Library/Fonts/RobotoCondensed-Regular.ttf",
                                bold = "/Users/zhangt8/Library/Fonts/RobotoCondensed-Bold.ttf",
                                italic = "/Users/zhangt8/Library/Fonts/RobotoCondensed-Italic.ttf",
                                bolditalic = "/Users/zhangt8/Library/Fonts/RobotoCondensed-BlackItalic.ttf",
                                dpi = 600) {
  font_files <- c(regular, bold, italic, bolditalic)
  if (!requireNamespace("showtext", quietly = TRUE) ||
      !requireNamespace("sysfonts", quietly = TRUE) ||
      !all(file.exists(font_files))) {
    options(clonetree.font_family = "")
    return(invisible(""))
  }
  sysfonts::font_add(
    family = family,
    regular = regular,
    bold = bold,
    italic = italic,
    bolditalic = bolditalic
  )
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = dpi)
  options(clonetree.font_family = family)
  options(clonetree.dpi = dpi)
  invisible(family)
}

ct_font_family <- function() {
  family <- getOption("clonetree.font_family", NULL)
  if (is.null(family)) {
    family <- use_clonetree_fonts()
  }
  family %||% ""
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "a", "alpha", "angle", "b", "ccf", "cluster_a", "cluster_b", "cluster_id",
    "fill_id", "hjust", "label", "label_x", "label_y", "mark", "outline",
    "outline_width", "relation", "sample_label", "row_id", "status", "x",
    "x_label", "x0", "x1", "xend", "xmax", "xmin", "y", "y0", "y1",
    "yend", "ymax", "ymin"
  ))
}

ct_arg_match <- function(x, choices) {
  x <- match.arg(x, choices)
  x
}

ct_first_non_missing <- function(x, default = NA) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    default
  } else {
    x[[1]]
  }
}

ct_is_numeric_id <- function(x) {
  all(grepl("^[0-9]+$", unique(as.character(x))))
}

ct_apply_cluster_prefix <- function(x, prefix = "C") {
  x <- as.character(x)
  if (is.null(prefix) || identical(prefix, "")) {
    return(x)
  }
  if (ct_is_numeric_id(x)) {
    paste0(prefix, x)
  } else {
    x
  }
}

ct_normalize_requested_clusters <- function(values, meta, prefix = "C") {
  if (is.null(values) || length(values) == 0) {
    return(character())
  }
  values <- unique(as.character(values))
  out <- character()
  for (value in values) {
    if (value %in% meta$cluster_id) {
      out <- c(out, value)
    } else if (value %in% meta$original_cluster_id) {
      out <- c(out, meta$cluster_id[match(value, meta$original_cluster_id)])
    } else if (!is.null(prefix) && startsWith(value, prefix)) {
      stripped <- sub(paste0("^", prefix), "", value)
      if (stripped %in% meta$original_cluster_id) {
        out <- c(out, meta$cluster_id[match(stripped, meta$original_cluster_id)])
      }
    } else if (!is.null(prefix)) {
      prefixed <- paste0(prefix, value)
      if (prefixed %in% meta$cluster_id) {
        out <- c(out, prefixed)
      }
    }
  }
  unique(out)
}

ct_strip_subject <- function(sample_ids, subject) {
  sample_ids <- as.character(sample_ids)
  if (is.null(subject) || length(subject) != 1 || is.na(subject)) {
    return(sample_ids)
  }
  prefix1 <- paste0(subject, "-")
  prefix2 <- paste0(subject, "_")
  out <- sample_ids
  out[startsWith(out, prefix1)] <- substring(out[startsWith(out, prefix1)], nchar(prefix1) + 1L)
  out[startsWith(out, prefix2)] <- substring(out[startsWith(out, prefix2)], nchar(prefix2) + 1L)
  out
}

ct_scale_ccf <- function(ccf, ccf_scale = c("auto", "fraction", "percent")) {
  ccf_scale <- match.arg(ccf_scale)
  ccf <- as.numeric(ccf)
  if (ccf_scale == "percent" || (ccf_scale == "auto" && suppressWarnings(max(ccf, na.rm = TRUE)) > 1.5)) {
    ccf <- ccf / 100
  }
  pmin(pmax(ccf, 0), 1)
}

ct_standardize_input <- function(data,
                                 subject = NULL,
                                 sample_col = "sampleID",
                                 cluster_col = "cluster.no",
                                 ccf_col = "CCF",
                                 mutation_col = "no.of.mutations",
                                 subject_col = "subject",
                                 clone_col = "Clone",
                                 trunk_value = "Y",
                                 cluster_prefix = "C",
                                 ccf_scale = c("auto", "fraction", "percent"),
                                 snap_ccf = TRUE,
                                 ccf_floor = 0.05,
                                 ccf_ceiling = 0.95,
                                 trunk_threshold = 0.95,
                                 summarize_ccf = c("max", "mean")) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  needed <- c(sample_col, cluster_col, ccf_col)
  missing <- setdiff(needed, names(data))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!is.null(subject_col) && subject_col %in% names(data)) {
    subjects <- unique(as.character(data[[subject_col]]))
    if (is.null(subject)) {
      if (length(subjects) > 1) {
        stop("`data` contains multiple subjects; provide `subject`.", call. = FALSE)
      }
      subject <- subjects[[1]]
    } else {
      data <- data[as.character(data[[subject_col]]) %in% as.character(subject), , drop = FALSE]
      if (nrow(data) == 0) {
        stop("No rows remain after filtering to `subject`.", call. = FALSE)
      }
      subject <- as.character(subject)[[1]]
    }
  }

  summarize_ccf <- match.arg(summarize_ccf)
  df <- data
  df$.sample_id <- as.character(df[[sample_col]])
  df$.sample_label <- ct_strip_subject(df$.sample_id, subject)
  df$.original_cluster_id <- as.character(df[[cluster_col]])
  df$.cluster_id <- ct_apply_cluster_prefix(df$.original_cluster_id, cluster_prefix)
  df$.ccf <- ct_scale_ccf(df[[ccf_col]], ccf_scale)
  df$.ccf[is.na(df$.ccf)] <- 0
  if (snap_ccf) {
    df$.ccf[df$.ccf <= ccf_floor] <- 0
    df$.ccf[df$.ccf >= ccf_ceiling] <- 1
  }

  cluster_order <- unique(df$.cluster_id)
  sample_order <- unique(df$.sample_id)
  sample_labels <- df$.sample_label[match(sample_order, df$.sample_id)]
  ccf_fun <- if (summarize_ccf == "max") max else mean
  agg <- stats::aggregate(
    .ccf ~ .cluster_id + .sample_id,
    data = df,
    FUN = ccf_fun
  )
  ccf <- matrix(
    0,
    nrow = length(cluster_order),
    ncol = length(sample_order),
    dimnames = list(cluster_order, sample_order)
  )
  ccf[cbind(match(agg$.cluster_id, cluster_order), match(agg$.sample_id, sample_order))] <- agg$.ccf

  if (mutation_col %in% names(df)) {
    mutations <- stats::aggregate(
      df[[mutation_col]],
      by = list(.cluster_id = df$.cluster_id),
      FUN = function(x) ct_first_non_missing(as.numeric(x), 1)
    )
    names(mutations)[2] <- "mutations"
  } else {
    mutations <- data.frame(.cluster_id = cluster_order, mutations = 1)
  }

  original <- stats::aggregate(
    df$.original_cluster_id,
    by = list(.cluster_id = df$.cluster_id),
    FUN = function(x) ct_first_non_missing(as.character(x), NA_character_)
  )
  names(original)[2] <- "original_cluster_id"

  if (!is.null(clone_col) && clone_col %in% names(df)) {
    trunk_df <- stats::aggregate(
      as.character(df[[clone_col]]) == trunk_value,
      by = list(.cluster_id = df$.cluster_id),
      FUN = any
    )
    names(trunk_df)[2] <- "is_trunk"
  } else {
    trunk_df <- data.frame(.cluster_id = cluster_order, is_trunk = rowMeans(ccf >= trunk_threshold) == 1)
  }

  meta <- merge(original, mutations, by = ".cluster_id", all = TRUE, sort = FALSE)
  meta <- merge(meta, trunk_df, by = ".cluster_id", all = TRUE, sort = FALSE)
  meta$cluster_id <- meta$.cluster_id
  meta$.cluster_id <- NULL
  meta$mean_ccf <- rowMeans(ccf[meta$cluster_id, , drop = FALSE])
  meta$max_ccf <- apply(ccf[meta$cluster_id, , drop = FALSE], 1, max)
  meta$is_trunk[is.na(meta$is_trunk)] <- FALSE
  meta$mutations[is.na(meta$mutations)] <- 1

  ord <- order(!meta$is_trunk, -meta$mean_ccf, -meta$mutations, meta$cluster_id)
  meta <- meta[ord, , drop = FALSE]
  rownames(meta) <- NULL
  ccf <- ccf[meta$cluster_id, , drop = FALSE]

  long <- data.frame(
    cluster_id = rep(rownames(ccf), times = ncol(ccf)),
    sample_id = rep(colnames(ccf), each = nrow(ccf)),
    sample_label = rep(sample_labels, each = nrow(ccf)),
    ccf = as.vector(ccf),
    stringsAsFactors = FALSE
  )
  long$original_cluster_id <- meta$original_cluster_id[match(long$cluster_id, meta$cluster_id)]
  long$mutations <- meta$mutations[match(long$cluster_id, meta$cluster_id)]
  long$is_trunk <- meta$is_trunk[match(long$cluster_id, meta$cluster_id)]

  structure(
    list(
      subject = subject,
      ccf = ccf,
      meta = meta,
      long = long,
      samples = data.frame(
        sample_id = sample_order,
        sample_label = sample_labels,
        stringsAsFactors = FALSE
      ),
      columns = list(
        sample_col = sample_col,
        cluster_col = cluster_col,
        ccf_col = ccf_col,
        mutation_col = mutation_col,
        subject_col = subject_col,
        clone_col = clone_col
      ),
      cluster_prefix = cluster_prefix
    ),
    class = "clonetree_input"
  )
}

ct_get_input <- function(x, ...) {
  if (inherits(x, "clonetree_input")) {
    x
  } else if (inherits(x, "clonetree_qc")) {
    x$input
  } else if (inherits(x, "clonetree")) {
    x$input
  } else {
    ct_standardize_input(x, ...)
  }
}

ct_save_plot <- function(plot, filename, width = NULL, height = NULL) {
  if (!is.null(filename)) {
    ct_font_family()
    ggplot2::ggsave(
      filename,
      plot = plot,
      width = width,
      height = height,
      dpi = getOption("clonetree.dpi", 600),
      bg = "white"
    )
    return(invisible(filename))
  }
  plot
}
