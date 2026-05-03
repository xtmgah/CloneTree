ct_cluster_colors <- function(clusters, trunk_clusters = character(), palette = clonetree_palette()) {
  clusters <- as.character(clusters)
  colors <- rep(palette, length.out = length(clusters) + 1L)
  values <- stats::setNames(colors[seq_along(clusters)], clusters)
  if (length(trunk_clusters) > 0) {
    values[intersect(trunk_clusters, names(values))] <- palette[[1]]
  }
  nontrunk <- setdiff(clusters, trunk_clusters)
  if (length(nontrunk) > 0) {
    values[nontrunk] <- rep(palette[-1], length.out = length(nontrunk))
  }
  values
}

ct_pretty_snvs <- function(x) {
  if (!is.finite(x) || x <= 0) {
    return(1)
  }
  exponent <- floor(log10(x))
  base <- x / 10^exponent
  nice <- c(1, 2, 5, 10)[which.min(abs(c(1, 2, 5, 10) - base))]
  nice * 10^exponent
}

ct_guess_event_label_col <- function(events) {
  candidates <- c("info", "event", "label", "driver", "alteration", "annotation")
  hit <- intersect(candidates, names(events))
  if (length(hit) == 0) {
    stop("Could not find an event label column. Provide `event_label_col`.", call. = FALSE)
  }
  hit[[1]]
}

ct_clean_event_label <- function(label) {
  lines <- unlist(strsplit(as.character(label), "\n", fixed = TRUE))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!tolower(lines) %in% c("branch event", "truncal driver")]
  paste(lines, collapse = "\n")
}

ct_prepare_events <- function(events,
                              tree,
                              layout,
                              event_cluster_col = "cluster.no",
                              event_label_col = NULL,
                              event_subject_col = "subject") {
  if (is.null(events)) {
    return(data.frame())
  }
  if (!is.data.frame(events)) {
    stop("`events` must be a data.frame.", call. = FALSE)
  }
  if (!event_cluster_col %in% names(events)) {
    alt <- intersect(c("cluster_id", "cluster", "cluster.no"), names(events))
    if (length(alt) == 0) {
      stop("Could not find an event cluster column. Provide `event_cluster_col`.", call. = FALSE)
    }
    event_cluster_col <- alt[[1]]
  }
  if (is.null(event_label_col)) {
    event_label_col <- ct_guess_event_label_col(events)
  }
  if (!event_label_col %in% names(events)) {
    stop("`event_label_col` is not present in `events`.", call. = FALSE)
  }
  if (!is.null(event_subject_col) && event_subject_col %in% names(events) && !is.null(tree$input$subject)) {
    events <- events[as.character(events[[event_subject_col]]) == tree$input$subject, , drop = FALSE]
  }
  if (nrow(events) == 0) {
    return(data.frame())
  }
  meta <- tree$input$meta
  event_cluster <- vapply(events[[event_cluster_col]], function(value) {
    normalized <- ct_normalize_requested_clusters(value, meta, tree$input$cluster_prefix)
    if (length(normalized) == 0) NA_character_ else normalized[[1]]
  }, character(1))
  keep <- event_cluster %in% tree$retained & !is.na(events[[event_label_col]]) & nzchar(as.character(events[[event_label_col]]))
  if (!any(keep)) {
    return(data.frame())
  }
  tmp <- data.frame(
    cluster_id = event_cluster[keep],
    label = vapply(events[[event_label_col]][keep], ct_clean_event_label, character(1)),
    stringsAsFactors = FALSE
  )
  tmp <- tmp[nzchar(tmp$label), , drop = FALSE]
  if (nrow(tmp) == 0) {
    return(data.frame())
  }
  tmp <- stats::aggregate(label ~ cluster_id, tmp, function(x) paste(unique(x), collapse = "\n"))
  edge_pos <- layout$edges[match(tmp$cluster_id, layout$edges$cluster_id), , drop = FALSE]
  tmp$x <- edge_pos$x + (edge_pos$xend - edge_pos$x) * 0.52
  tmp$y <- edge_pos$y + (edge_pos$yend - edge_pos$y) * 0.52
  dx <- edge_pos$xend - edge_pos$x
  dy <- edge_pos$yend - edge_pos$y
  norm <- sqrt(dx^2 + dy^2)
  norm[norm == 0] <- 1
  tmp$perp_x <- -dy / norm
  tmp$perp_y <- dx / norm
  tmp$nudge_x <- ifelse(edge_pos$xend >= edge_pos$x, 0.34, -0.34)
  tmp$nudge_y <- ifelse(edge_pos$yend <= edge_pos$y, 0.09, -0.09)
  tmp$is_cna <- grepl("ASCNA|BCNA|DLOH|NLOH|HOMD|amp|del|gain|loss", tmp$label, ignore.case = TRUE)
  tmp
}

ct_oval_cluster_order <- function(input, clusters, presence_threshold = 0.02) {
  groups <- ct_oval_cluster_groups(input, clusters, presence_threshold)
  c(groups$trunk, groups$branches, groups$leaves)
}

ct_oval_cluster_groups <- function(input, clusters, presence_threshold = 0.02) {
  clusters <- intersect(input$meta$cluster_id, clusters)
  meta <- input$meta[match(clusters, input$meta$cluster_id), , drop = FALSE]
  trunk <- clusters[meta$is_trunk]
  subclones <- setdiff(clusters, trunk)
  ccf <- input$ccf[subclones, , drop = FALSE]
  present_count <- if (length(subclones) == 0) {
    numeric()
  } else {
    rowSums(ccf > presence_threshold)
  }
  branches <- subclones[present_count > 1]
  leaves <- subclones[present_count <= 1]
  list(trunk = trunk, branches = branches, leaves = leaves)
}

ct_oval_labels <- function(input, clusters, status = NULL) {
  meta <- input$meta[match(clusters, input$meta$cluster_id), , drop = FALSE]
  total <- sum(meta$mutations, na.rm = TRUE)
  if (!is.finite(total) || total == 0) {
    total <- 1
  }
  labels <- paste0(
    meta$cluster_id,
    "\n",
    meta$mutations,
    "\n",
    sprintf("%.1f%%", 100 * meta$mutations / total)
  )
  if (!is.null(status)) {
    labels[status[meta$cluster_id] == "excluded"] <- paste0(labels[status[meta$cluster_id] == "excluded"], "\nexcluded")
  }
  stats::setNames(labels, meta$cluster_id)
}

ct_oval_plot <- function(qc,
                         solution = 1,
                         show_values = TRUE,
                         clonesum = FALSE,
                         include_excluded = TRUE,
                         solution_position = c("subtitle", "caption"),
                         subclone_outline = FALSE,
                         summary_subclone_alpha = 0.5,
                         add_group_spacing = FALSE,
                         text_size = 11,
                         palette = clonetree_palette(),
                         oval_width = 6,
                         oval_height = 10,
                         panel_fill = "#EBEBEA",
                         filename = NULL,
                         width = NULL,
                         height = NULL) {
  if (solution < 1 || solution > length(qc$solutions)) {
    stop("`solution` must be between 1 and ", length(qc$solutions), ".", call. = FALSE)
  }
  solution_position <- match.arg(solution_position)
  input <- qc$input
  retained <- qc$solutions[[solution]]$retained
  solution_excluded <- qc$solutions[[solution]]$excluded
  active <- if (include_excluded) qc$active_clusters else retained
  excluded <- setdiff(active, retained)
  cluster_groups <- ct_oval_cluster_groups(input, active, qc$presence_threshold)
  cluster_order <- c(cluster_groups$trunk, cluster_groups$branches, cluster_groups$leaves)
  if (length(cluster_order) == 0) {
    stop("No clusters are available to plot.", call. = FALSE)
  }

  df <- input$long[input$long$cluster_id %in% cluster_order, , drop = FALSE]
  df$is_trunk <- input$meta$is_trunk[match(df$cluster_id, input$meta$cluster_id)]
  df$status <- ifelse(df$cluster_id %in% retained, "retained", "excluded")
  status <- stats::setNames(ifelse(cluster_order %in% retained, "retained", "excluded"), cluster_order)
  df$fill_id <- ifelse(df$status == "retained", df$cluster_id, "excluded")
  df$x_offset <- 0
  df$y_offset <- 0
  df$angle <- 0
  df$row_id <- df$cluster_id
  df$alpha <- 1
  df$outline <- ifelse(subclone_outline & !df$is_trunk & df$status == "retained", "subclone", df$status)
  df$outline_width <- ifelse(df$outline == "subclone", 0.25, ifelse(df$status == "excluded", 0.45, 0.1))

  row_levels <- cluster_order
  if (clonesum) {
    clone_sum_id <- "Clone summary"
    set.seed(sum(utf8ToInt(input$subject %||% "CloneTree")) %% .Machine$integer.max)
    summary_df <- df[df$status == "retained", , drop = FALSE]
    if (nrow(summary_df) > 0) {
      summary_df$row_id <- clone_sum_id
      nonzero <- summary_df$ccf > 0
      summary_df$x_offset[nonzero] <- (1 - summary_df$ccf[nonzero]) * stats::rnorm(sum(nonzero), 0, 0.08)
      summary_df$y_offset[nonzero] <- (1 - summary_df$ccf[nonzero]) * stats::rnorm(sum(nonzero), 0, 0.20)
      summary_df$angle[nonzero] <- (1 - summary_df$ccf[nonzero]) * stats::rnorm(sum(nonzero), 0, 6)
      summary_df$alpha[!summary_df$is_trunk] <- summary_subclone_alpha
      df <- rbind(summary_df, df[df$status == "retained", , drop = FALSE])
      row_levels <- c(clone_sum_id, cluster_order[cluster_order %in% retained])
    }
  }

  df$row_id <- factor(df$row_id, levels = row_levels)
  df$sample_label <- factor(df$sample_label, levels = input$samples$sample_label)
  colors <- ct_cluster_colors(cluster_order, input$meta$cluster_id[input$meta$is_trunk], palette)
  colors <- c(colors, excluded = "#C9C9C9")
  row_labels <- ct_oval_labels(input, cluster_order, status)
  if (clonesum) {
    row_labels <- c(
      "Clone summary" = "Clone\nsummary",
      row_labels[names(row_labels) %in% row_levels]
    )
  }
  solution_label <- if (length(solution_excluded) == 0) {
    paste0("Solution ", solution, "/", length(qc$solutions), ": no clusters excluded")
  } else {
    paste0("Solution ", solution, "/", length(qc$solutions), ": excluded ", paste(solution_excluded, collapse = ", "))
  }

  font_family <- ct_font_family()
  sample_levels <- levels(df$sample_label)
  row_group <- stats::setNames(rep("branch", length(row_levels)), row_levels)
  if (clonesum) {
    row_group[clone_sum_id] <- "summary"
  }
  row_group[intersect(row_levels, cluster_groups$trunk)] <- "trunk"
  row_group[intersect(row_levels, cluster_groups$branches)] <- "branch"
  row_group[intersect(row_levels, cluster_groups$leaves)] <- "leaf"
  row_step <- 1.48
  group_gap <- if (add_group_spacing) 0.12 else 0
  row_y <- numeric(length(row_levels))
  cursor <- 0
  for (i in seq_along(row_levels)) {
    if (i > 1) {
      cursor <- cursor - row_step
      if (!identical(row_group[[row_levels[[i]]]], row_group[[row_levels[[i - 1L]]]])) {
        cursor <- cursor - group_gap
      }
    }
    row_y[[i]] <- cursor
  }
  row_y <- stats::setNames(row_y, row_levels)
  sample_x <- stats::setNames(seq_along(sample_levels), sample_levels)
  ellipse_a <- 0.40
  ellipse_b <- ellipse_a * oval_height / oval_width
  cell_half_w <- 0.48
  cell_half_h <- max(0.70, ellipse_b * 1.05)

  df$x0 <- unname(sample_x[as.character(df$sample_label)]) + df$x_offset
  df$y0 <- unname(row_y[as.character(df$row_id)]) + df$y_offset
  df$a <- ellipse_a * df$ccf
  df$b <- ellipse_b * df$ccf

  bg_df <- expand.grid(
    row_id = factor(row_levels, levels = row_levels),
    sample_label = factor(sample_levels, levels = sample_levels),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  bg_df$x <- unname(sample_x[as.character(bg_df$sample_label)])
  bg_df$y <- unname(row_y[as.character(bg_df$row_id)])
  bg_df$xmin <- bg_df$x - cell_half_w
  bg_df$xmax <- bg_df$x + cell_half_w
  bg_df$ymin <- bg_df$y - cell_half_h
  bg_df$ymax <- bg_df$y + cell_half_h

  label_df <- data.frame(
    row_id = factor(row_levels, levels = row_levels),
    label = unname(row_labels[row_levels]),
    x = length(sample_levels) + 0.68,
    y = unname(row_y[row_levels]),
    stringsAsFactors = FALSE
  )
  sample_df <- data.frame(
    sample_label = sample_levels,
    label = sample_levels,
    x = unname(sample_x[sample_levels]),
    y = max(row_y) + cell_half_h + 0.22,
    stringsAsFactors = FALSE
  )
  plot_title <- if (is.null(input$subject)) {
    if (show_values) "CloneTree CCF QC" else "CloneTree"
  } else if (show_values) {
    paste0(input$subject, " CCF QC")
  } else {
    input$subject
  }
  has_top_solution <- solution_position == "subtitle"
  title_df <- data.frame(
    x = mean(unname(sample_x[sample_levels])),
    y = max(sample_df$y) + if (has_top_solution) 0.78 else 0.58,
    label = plot_title,
    stringsAsFactors = FALSE
  )
  subtitle_df <- data.frame(
    x = title_df$x,
    y = max(sample_df$y) + 0.42,
    label = solution_label,
    stringsAsFactors = FALSE
  )
  caption_df <- data.frame(
    x = title_df$x,
    y = min(row_y) - cell_half_h - 0.18,
    label = solution_label,
    stringsAsFactors = FALSE
  )
  y_min <- if (solution_position == "caption") caption_df$y - 0.12 else min(row_y) - cell_half_h - 0.20
  y_limits <- c(y_min, max(title_df$y) + 0.18)
  x_limits <- c(0.48, length(sample_levels) + 1.78)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bg_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = panel_fill,
      color = "white",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggforce::geom_ellipse(
      data = df,
      ggplot2::aes(
        x0 = x0,
        y0 = y0,
        a = a,
        b = b,
        angle = angle,
        fill = fill_id,
        color = outline,
        alpha = alpha,
        linewidth = outline_width
      ),
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_color_manual(
      values = c(retained = "transparent", excluded = "#B2182B", subclone = "black"),
      guide = "none"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::geom_text(
      data = title_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = text_size * 0.45,
      family = font_family,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = sample_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = text_size * 0.36,
      family = font_family,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0,
      lineheight = 0.9,
      size = text_size * 0.30,
      family = font_family,
      inherit.aes = FALSE
    ) +
    ggplot2::coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
    ggplot2::theme_void(base_size = text_size, base_family = font_family) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = grid::unit(c(6, 18, 8, 8), "pt")
    )

  if (has_top_solution) {
    p <- p + ggplot2::geom_text(
      data = subtitle_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = text_size * 0.34,
      family = font_family,
      inherit.aes = FALSE
    )
  }

  if (solution_position == "caption") {
    p <- p + ggplot2::geom_text(
      data = caption_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = text_size * 0.34,
      family = font_family,
      inherit.aes = FALSE
    )
  }

  if (show_values) {
    value_df <- df[!clonesum | df$row_id != "Clone summary", , drop = FALSE]
    p <- p + ggplot2::geom_text(
      data = value_df,
      ggplot2::aes(x = x0, y = y0, label = sprintf("%.2f", ccf)),
      size = text_size * 0.26,
      family = font_family,
      inherit.aes = FALSE
    )
  }

  ct_save_plot(p, filename, width, height)
}

#' Plot a publication-style clonal phylogenetic tree
#'
#' @param x A `clonetree`, `clonetree_qc`, `clonetree_input`, or raw clustering
#'   data frame.
#' @param solution Which QC solution to draw when `x` is not already a tree.
#' @param events Optional data frame of driver mutations or copy-number events
#'   to mark on branches.
#' @param event_cluster_col,event_label_col,event_subject_col Column names in
#'   `events`. If `event_label_col` is `NULL`, CloneTree looks for `info`,
#'   `event`, `label`, `driver`, `alteration`, or `annotation`.
#' @param title Plot title. Defaults to the subject ID when available.
#' @param palette Cluster color palette.
#' @param branch_width,node_size Line width and node size.
#' @param branch_length `"linear"` or `"sqrt"` visual branch-length transform.
#' @param min_branch_fraction Minimum visual branch length as a fraction of the
#'   longest branch.
#' @param angle_width Fan width below the trunk. The default keeps first-level
#'   branches within about 35 degrees of the vertical trunk.
#' @param show_scale,label_samples,label_events Logical toggles.
#' @param filename Optional output path. When supplied, the plot is saved and
#'   returned invisibly.
#' @param width,height Optional output dimensions for `filename`.
#' @param ... Passed to [infer_clone_tree()] when needed.
#'
#' @return A `ggplot` object, or invisibly the filename if saved.
#' @export
plot_clone_tree <- function(x,
                            solution = 1,
                            events = NULL,
                            event_cluster_col = "cluster.no",
                            event_label_col = NULL,
                            event_subject_col = "subject",
                            title = NULL,
                            palette = clonetree_palette(),
                            branch_width = 3.2,
                            node_size = 4.4,
                            branch_length = c("linear", "sqrt"),
                            min_branch_fraction = 0.08,
                            angle_width = 70,
                            show_scale = TRUE,
                            label_samples = TRUE,
                            label_events = TRUE,
                            filename = NULL,
                            width = NULL,
                            height = NULL,
                            ...) {
  tree <- if (inherits(x, "clonetree")) {
    x
  } else {
    infer_clone_tree(x, solution = solution, ...)
  }
  branch_length <- match.arg(branch_length)
  layout <- layout_clone_tree(
    tree,
    branch_length = branch_length,
    min_branch_fraction = min_branch_fraction,
    angle_width = angle_width
  )
  colors <- ct_cluster_colors(tree$retained, tree$trunk_clusters, palette)
  nodes <- layout$nodes
  edges <- layout$edges
  title <- title %||% tree$input$subject %||% ""
  font_family <- ct_font_family()

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = cluster_id),
      linewidth = branch_width,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y),
      color = "#D9D9D9",
      fill = "#D9D9D9",
      size = node_size
    ) +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void(base_size = 11, base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = 13, family = font_family),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = grid::unit(c(8, 24, 8, 24), "pt")
    ) +
    ggplot2::labs(title = title)

  gl <- nodes[nodes$node_id == "GL", , drop = FALSE]
  if (nrow(gl) == 1) {
    x_span_gl <- diff(range(nodes$x, na.rm = TRUE))
    if (!is.finite(x_span_gl) || x_span_gl == 0) x_span_gl <- 1
    x_span_gl <- max(x_span_gl, 1)
    p <- p + ggplot2::geom_text(
      data = gl,
      ggplot2::aes(x = x - 0.045 * x_span_gl, y = y, label = "GL"),
      hjust = 1,
      vjust = 0.5,
      size = 3.6,
      family = font_family
    )
  }

  if (label_samples) {
    leaf <- nodes[nodes$is_terminal & nzchar(nodes$sample_label), , drop = FALSE]
    if (nrow(leaf) > 0) {
      span <- diff(range(nodes$x, na.rm = TRUE))
      if (!is.finite(span) || span == 0) span <- 1
      span <- max(span, 1)
      leaf$label_x <- leaf$x + ifelse(leaf$x >= 0, 0.09, -0.09) * span
      leaf$label_y <- leaf$y - 0.01
      leaf$hjust <- ifelse(leaf$x >= 0, 0, 1)
      p <- p + ggplot2::geom_text(
        data = leaf,
        ggplot2::aes(x = label_x, y = label_y, label = sample_label, hjust = hjust),
        size = 3.4,
        fontface = "bold",
        family = font_family,
        inherit.aes = FALSE
      )
    }
  }

  event_df <- ct_prepare_events(
    events = events,
    tree = tree,
    layout = layout,
    event_cluster_col = event_cluster_col,
    event_label_col = event_label_col,
    event_subject_col = event_subject_col
  )
  if (label_events && nrow(event_df) > 0) {
    span <- diff(range(nodes$x, event_df$x, na.rm = TRUE))
    if (!is.finite(span) || span == 0) span <- 1
    span <- max(span, 1)
    p <- p + ggrepel::geom_text_repel(
      data = event_df,
      ggplot2::aes(x = x, y = y, label = label, color = cluster_id),
      nudge_x = event_df$nudge_x * span,
      nudge_y = event_df$nudge_y,
      min.segment.length = 0,
      segment.size = 0.25,
      box.padding = 0.35,
      point.padding = 0.25,
      force = 8,
      force_pull = 0,
      max.overlaps = Inf,
      size = 3.1,
      lineheight = 0.95,
      seed = 1,
      family = font_family,
      inherit.aes = FALSE
    )
  }

  if (show_scale && nrow(edges) > 0) {
    max_mut <- max(edges$mutations, na.rm = TRUE)
    scale_snvs <- ct_pretty_snvs(max_mut / 4)
    if (branch_length == "sqrt") {
      scale_len <- sqrt(scale_snvs) / sqrt(max_mut)
    } else {
      scale_len <- scale_snvs / max_mut
    }
    x_span <- diff(range(nodes$x, na.rm = TRUE))
    if (!is.finite(x_span) || x_span == 0) x_span <- 1
    x_span <- max(x_span, 1)
    scale_x <- max(nodes$x, na.rm = TRUE) + 0.32 * x_span
    scale_y <- max(nodes$y, na.rm = TRUE) - scale_len
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(x = scale_x, xend = scale_x, y = scale_y, yend = scale_y + scale_len),
        linewidth = 0.9,
        color = "black",
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = scale_x, y = scale_y + scale_len / 2, label = paste0(scale_snvs, " SNVs")),
        hjust = 1.15,
        size = 3,
        family = font_family,
        inherit.aes = FALSE
      )
  }

  ct_save_plot(p, filename, width, height)
}

#' Plot the QC pairwise rule matrix
#'
#' @param qc A `clonetree_qc` object.
#' @param ccf_label Which cluster-level CCF summary to show under each cluster
#'   name on both axes.
#' @param filename,width,height Optional save arguments.
#'
#' @return A `ggplot` object.
#' @export
plot_clone_qc <- function(qc,
                          ccf_label = c("max", "mean", "none"),
                          filename = NULL,
                          width = NULL,
                          height = NULL) {
  if (!inherits(qc, "clonetree_qc")) {
    stop("`qc` must be a clonetree_qc object.", call. = FALSE)
  }
  ccf_label <- match.arg(ccf_label)
  pairwise <- qc$pairwise
  if (nrow(pairwise) == 0) {
    stop("No pairwise rules to plot.", call. = FALSE)
  }
  meta <- qc$input$meta[match(qc$active_clusters, qc$input$meta$cluster_id), , drop = FALSE]
  axis_labels <- qc$active_clusters
  if (ccf_label != "none") {
    ccf_values <- if (ccf_label == "max") meta$max_ccf else meta$mean_ccf
    axis_labels <- paste0(qc$active_clusters, "\n", sprintf("%.2f", ccf_values))
  }
  axis_labels <- stats::setNames(axis_labels, qc$active_clusters)
  levels <- c("unconstrained", "nested_allowed", "branch_required", "nested_required", "incompatible")
  labels <- c("No rule", "Nested allowed", "Branch required", "Nested required", "Exclude one")
  pairwise$relation <- factor(pairwise$relation, levels = levels, labels = labels)
  pairwise$cluster_a <- factor(pairwise$cluster_a, levels = qc$active_clusters)
  pairwise$cluster_b <- factor(pairwise$cluster_b, levels = rev(qc$active_clusters))
  pairwise$mark <- ifelse(pairwise$incompatible, "PH+X", ifelse(pairwise$pigeonhole, "PH", ifelse(pairwise$crossing, "X", "")))

  font_family <- ct_font_family()
  p <- ggplot2::ggplot(pairwise, ggplot2::aes(x = cluster_a, y = cluster_b, fill = relation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = mark), size = 3, family = font_family) +
    ggplot2::scale_fill_manual(
      values = c(
        "No rule" = "#EFEFEF",
        "Nested allowed" = "#D9E8C8",
        "Branch required" = "#C9DCEB",
        "Nested required" = "#F2D49B",
        "Exclude one" = "#D95F5F"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::scale_y_discrete(labels = axis_labels) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 11, base_family = font_family) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, lineheight = 0.9, family = font_family),
      axis.text.y = ggplot2::element_text(lineheight = 0.9, family = font_family),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(family = font_family),
      plot.subtitle = ggplot2::element_text(family = font_family)
    ) +
    ggplot2::labs(
      title = if (is.null(qc$input$subject)) "CloneTree QC rules" else paste0(qc$input$subject, " QC rules"),
      subtitle = if (ccf_label == "none") {
        "PH = pigeonhole; X = crossing"
      } else {
        paste0("PH = pigeonhole; X = crossing; axis value = ", ccf_label, " CCF")
      }
    )

  ct_save_plot(p, filename, width, height)
}

#' Plot cluster CCF ovals for QC review
#'
#' @param qc A `clonetree_qc` object.
#' @param solution Which solution to highlight.
#' @param show_values If `TRUE`, print rounded CCF values inside ovals.
#' @param clonesum If `TRUE`, add a top clone-summary row with retained
#'   clusters overlaid within each region, matching the original draft plot.
#' @param include_excluded If `TRUE`, show excluded clusters in grey with a red
#'   outline. Use `FALSE` for the final retained-cluster plot.
#' @param solution_position Where to print the solution summary: `"subtitle"`
#'   or `"caption"`.
#' @param subclone_outline If `TRUE`, draw a thin black outline around retained
#'   subclone ovals.
#' @param summary_subclone_alpha Fill alpha for subclone ovals in the clone
#'   summary row.
#' @param add_group_spacing If `TRUE`, increase vertical white space between
#'   clone-summary, trunk, branch, and leaf sections.
#' @param text_size Base text size for oval plots.
#' @param palette Cluster colors.
#' @param oval_width,oval_height Ellipse radii multipliers, as in the original
#'   draft plotting function.
#' @param filename,width,height Optional save arguments.
#'
#' @return A `ggplot` object.
#' @export
plot_clone_ccf <- function(qc,
                           solution = 1,
                           show_values = TRUE,
                           clonesum = FALSE,
                           include_excluded = TRUE,
                           solution_position = c("subtitle", "caption"),
                           subclone_outline = FALSE,
                           summary_subclone_alpha = 0.5,
                           add_group_spacing = FALSE,
                           text_size = 11,
                           palette = clonetree_palette(),
                           oval_width = 6,
                           oval_height = 10,
                           filename = NULL,
                           width = NULL,
                           height = NULL) {
  if (!inherits(qc, "clonetree_qc")) {
    stop("`qc` must be a clonetree_qc object.", call. = FALSE)
  }
  ct_oval_plot(
    qc = qc,
    solution = solution,
    show_values = show_values,
    clonesum = clonesum,
    include_excluded = include_excluded,
    solution_position = solution_position,
    subclone_outline = subclone_outline,
    summary_subclone_alpha = summary_subclone_alpha,
    add_group_spacing = add_group_spacing,
    text_size = text_size,
    palette = palette,
    oval_width = oval_width,
    oval_height = oval_height,
    filename = filename,
    width = width,
    height = height
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}
