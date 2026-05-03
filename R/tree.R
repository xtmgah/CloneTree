ct_crosses <- function(a, b, tolerance) {
  diff <- a - b
  any(diff > tolerance) && any(diff < -tolerance)
}

ct_can_parent <- function(parent_ccf, child_ccf, tolerance) {
  all(parent_ccf + tolerance >= child_ccf) && any(parent_ccf > child_ccf + tolerance)
}

ct_is_ancestor <- function(parent_map, ancestor, node) {
  current <- parent_map[[node]]
  while (!is.null(current) && !is.na(current) && current != "GL") {
    if (identical(current, ancestor)) {
      return(TRUE)
    }
    current <- parent_map[[current]]
  }
  FALSE
}

ct_parent_depth <- function(parent_map, node) {
  depth <- 0L
  current <- parent_map[[node]]
  while (!is.null(current) && !is.na(current)) {
    depth <- depth + 1L
    if (current == "GL") {
      break
    }
    current <- parent_map[[current]]
  }
  depth
}

ct_order_edges <- function(edges) {
  if (nrow(edges) == 0) {
    return(edges)
  }
  depths <- vapply(edges$child, function(x) ct_parent_depth(stats::setNames(edges$parent, edges$child), x), integer(1))
  edges[order(depths, !edges$is_trunk, -edges$mutations, edges$child), , drop = FALSE]
}

ct_infer_parent_map <- function(input,
                                retained,
                                pairwise,
                                parent_tolerance = 0.03,
                                presence_threshold = 0.02) {
  meta <- input$meta
  ccf <- input$ccf[retained, , drop = FALSE]
  trunks <- retained[retained %in% meta$cluster_id[meta$is_trunk]]
  nontrunks <- setdiff(retained, trunks)
  parent <- stats::setNames(rep(NA_character_, length(retained)), retained)

  if (length(trunks) > 0) {
    parent[trunks[1]] <- "GL"
    if (length(trunks) > 1) {
      for (i in 2:length(trunks)) {
        parent[trunks[i]] <- trunks[i - 1L]
      }
    }
    default_parent <- trunks[length(trunks)]
  } else {
    default_parent <- "GL"
  }

  mean_ccf <- stats::setNames(meta$mean_ccf[match(retained, meta$cluster_id)], retained)
  for (child in nontrunks) {
    child_ccf <- ccf[child, ]
    candidates <- setdiff(retained, child)
    candidates <- candidates[vapply(candidates, function(candidate) {
      ct_can_parent(ccf[candidate, ], child_ccf, parent_tolerance)
    }, logical(1))]
    if (length(candidates) > 0) {
      score <- vapply(candidates, function(candidate) {
        diff <- ccf[candidate, ] - child_ccf
        present <- child_ccf > presence_threshold
        if (any(present)) {
          mean(diff[present])
        } else {
          mean(diff)
        }
      }, numeric(1))
      candidates <- candidates[order(score, abs(mean_ccf[candidates] - mean_ccf[child]), candidates)]
      parent[child] <- candidates[[1]]
    } else {
      parent[child] <- default_parent
    }
  }

  nested <- pairwise[pairwise$relation == "nested_required", , drop = FALSE]
  if (nrow(nested) > 0) {
    for (i in seq_len(nrow(nested))) {
      a <- nested$cluster_a[i]
      b <- nested$cluster_b[i]
      if (!all(c(a, b) %in% retained)) {
        next
      }
      if (nested$a_can_parent_b[i] && !ct_is_ancestor(parent, a, b)) {
        if (!ct_is_ancestor(parent, b, a)) {
          parent[b] <- a
        }
      } else if (nested$b_can_parent_a[i] && !ct_is_ancestor(parent, b, a)) {
        if (!ct_is_ancestor(parent, a, b)) {
          parent[a] <- b
        }
      }
    }
  }

  parent
}

#' Infer a clone tree from a QC solution
#'
#' @description
#' Builds a deterministic parent-child hierarchy from one retained-cluster
#' solution. Parent candidates must have CCF greater than or equal to the child
#' in every sample, within tolerance. The closest such ancestor is selected, so
#' compatible pigeonhole pairs are automatically nested under the trunk instead
#' of being drawn as disconnected bars.
#'
#' @param x A raw clustering data frame, `clonetree_input`, or `clonetree_qc`.
#' @param solution Which QC solution to use.
#' @param parent_tolerance Tolerance for accepting a parent CCF as greater than
#'   or equal to a child CCF.
#' @param presence_threshold CCF threshold used to label cluster/sample
#'   presence at terminal nodes.
#' @param ... Passed to [clone_qc()] when `x` is not already a QC result.
#'
#' @return A `clonetree` object.
#' @export
infer_clone_tree <- function(x,
                             solution = 1,
                             parent_tolerance = 0.03,
                             presence_threshold = 0.02,
                             ...) {
  qc <- if (inherits(x, "clonetree_qc")) {
    x
  } else {
    clone_qc(x, presence_threshold = presence_threshold, ...)
  }
  if (solution < 1 || solution > length(qc$solutions)) {
    stop("`solution` must be between 1 and ", length(qc$solutions), ".", call. = FALSE)
  }
  input <- qc$input
  retained <- qc$solutions[[solution]]$retained
  meta <- input$meta[match(retained, input$meta$cluster_id), , drop = FALSE]
  parent <- ct_infer_parent_map(
    input = input,
    retained = retained,
    pairwise = qc$pairwise,
    parent_tolerance = parent_tolerance,
    presence_threshold = presence_threshold
  )

  edges <- data.frame(
    parent = unname(parent),
    child = names(parent),
    cluster_id = names(parent),
    original_cluster_id = meta$original_cluster_id[match(names(parent), meta$cluster_id)],
    mutations = meta$mutations[match(names(parent), meta$cluster_id)],
    mean_ccf = meta$mean_ccf[match(names(parent), meta$cluster_id)],
    max_ccf = meta$max_ccf[match(names(parent), meta$cluster_id)],
    is_trunk = meta$is_trunk[match(names(parent), meta$cluster_id)],
    stringsAsFactors = FALSE
  )
  edges <- ct_order_edges(edges)

  node_ids <- c("GL", retained)
  nodes <- data.frame(
    node_id = node_ids,
    cluster_id = c(NA_character_, retained),
    label = node_ids,
    is_germline = node_ids == "GL",
    is_trunk = c(FALSE, meta$is_trunk[match(retained, meta$cluster_id)]),
    stringsAsFactors = FALSE
  )

  ccf <- input$ccf[retained, , drop = FALSE]
  present_samples <- lapply(retained, function(cluster) {
    labs <- input$samples$sample_label[match(colnames(ccf), input$samples$sample_id)]
    labs[ccf[cluster, ] > presence_threshold]
  })
  names(present_samples) <- retained

  out <- list(
    input = input,
    qc = qc,
    solution = solution,
    retained = retained,
    excluded = qc$solutions[[solution]]$excluded,
    parent = parent,
    nodes = nodes,
    edges = edges,
    present_samples = present_samples,
    parent_tolerance = parent_tolerance,
    presence_threshold = presence_threshold,
    trunk_clusters = retained[retained %in% meta$cluster_id[meta$is_trunk]]
  )
  class(out) <- "clonetree"
  out
}

ct_scaled_lengths <- function(edges,
                              branch_length = c("linear", "sqrt"),
                              min_branch_fraction = 0.04) {
  branch_length <- match.arg(branch_length)
  raw <- edges$mutations
  raw[is.na(raw) | raw < 0] <- 0
  scaled <- if (branch_length == "sqrt") sqrt(raw) else raw
  if (max(scaled, na.rm = TRUE) == 0) {
    scaled <- rep(1, length(scaled))
  } else {
    scaled <- scaled / max(scaled, na.rm = TRUE)
  }
  pmax(scaled, min_branch_fraction)
}

ct_terminal_nodes <- function(edges, retained) {
  setdiff(retained, edges$parent[edges$parent %in% retained])
}

ct_sibling_angles <- function(kids,
                              is_terminal,
                              center_angle = 0,
                              width = 70,
                              max_abs_angle = 60) {
  n_kids <- length(kids)
  angles <- stats::setNames(rep(0, n_kids), kids)
  if (n_kids == 0) {
    return(angles)
  }
  if (n_kids == 1) {
    angles[[kids[[1]]]] <- if (is_terminal[[1]]) 0 else center_angle
    return(angles)
  }

  fan_width <- min(width, max(50, 28 * (n_kids - 1)))
  angle_center <- if (all(is_terminal)) 0 else center_angle
  angle_center <- max(
    min(angle_center, max_abs_angle - fan_width / 2),
    -max_abs_angle + fan_width / 2
  )
  slots <- seq(angle_center - fan_width / 2, angle_center + fan_width / 2, length.out = n_kids)

  nonterminal <- kids[!is_terminal]
  terminal <- kids[is_terminal]
  center_slots <- order(abs(slots - angle_center), slots)
  assigned <- integer(0)
  if (length(nonterminal) > 0) {
    keep <- center_slots[seq_len(length(nonterminal))]
    angles[nonterminal] <- slots[keep]
    assigned <- c(assigned, keep)
  }
  remaining_slots <- setdiff(seq_along(slots), assigned)
  if (length(terminal) > 0) {
    angles[terminal] <- slots[remaining_slots]
  }
  angles
}

#' Compute plotting coordinates for a clone tree
#'
#' @param tree A `clonetree` object.
#' @param branch_length `"linear"` or `"sqrt"` transformation before plotting.
#' @param min_branch_fraction Minimum visual branch length as a fraction of the
#'   longest plotted branch.
#' @param angle_width Fan width, in degrees, for branches below the trunk.
#'   The default keeps first-level branches close to the vertical trunk and
#'   narrows descendant branches at deeper levels.
#'
#' @return A list with `nodes` and `edges` data frames containing coordinates.
#' @export
layout_clone_tree <- function(tree,
                              branch_length = c("linear", "sqrt"),
                              min_branch_fraction = 0.08,
                              angle_width = 70) {
  if (!inherits(tree, "clonetree")) {
    stop("`tree` must be a clonetree object.", call. = FALSE)
  }
  branch_length <- match.arg(branch_length)
  edges <- tree$edges
  edges$.length <- ct_scaled_lengths(edges, branch_length, min_branch_fraction)

  children <- split(edges$child, edges$parent)
  edge_lookup <- stats::setNames(seq_len(nrow(edges)), edges$child)
  node_pos <- data.frame(
    node_id = "GL",
    x = 0,
    y = 0,
    stringsAsFactors = FALSE
  )
  edge_pos <- edges
  edge_pos$x <- edge_pos$y <- edge_pos$xend <- edge_pos$yend <- NA_real_

  trunk_chain <- tree$trunk_clusters
  trunk_lengths <- if (length(trunk_chain) > 0) {
    edges$.length[match(trunk_chain, edges$child)]
  } else {
    numeric()
  }
  current <- "GL"
  x <- 0
  y <- sum(trunk_lengths)
  node_pos$y[node_pos$node_id == "GL"] <- y

  if (length(trunk_chain) > 0) {
    for (trunk in trunk_chain) {
      idx <- edge_lookup[[trunk]]
      len <- edges$.length[idx]
      edge_pos$x[idx] <- x
      edge_pos$y[idx] <- y
      edge_pos$xend[idx] <- x
      edge_pos$yend[idx] <- y - len
      y <- y - len
      node_pos <- rbind(
        node_pos,
        data.frame(node_id = trunk, x = x, y = y, stringsAsFactors = FALSE)
      )
      current <- trunk
    }
  }

  place_children <- function(parent, px, py, center_angle = 0, width = angle_width, depth = 1) {
    kids <- children[[parent]]
    kids <- setdiff(kids, trunk_chain)
    if (length(kids) == 0) {
      return(invisible(NULL))
    }
    kid_meta <- tree$input$meta[match(kids, tree$input$meta$cluster_id), , drop = FALSE]
    kids <- kids[order(-kid_meta$max_ccf, -kid_meta$mutations, kids)]
    is_terminal <- !kids %in% names(children)
    angles <- ct_sibling_angles(kids, is_terminal, center_angle = center_angle, width = width)
    for (i in seq_along(kids)) {
      child <- kids[[i]]
      idx <- edge_lookup[[child]]
      len <- edges$.length[idx]
      theta <- angles[[child]] * pi / 180
      cx <- px + len * sin(theta)
      cy <- py - len * cos(theta)
      edge_pos$x[idx] <<- px
      edge_pos$y[idx] <<- py
      edge_pos$xend[idx] <<- cx
      edge_pos$yend[idx] <<- cy
      node_pos <<- rbind(
        node_pos,
        data.frame(node_id = child, x = cx, y = cy, stringsAsFactors = FALSE)
      )
      place_children(
        parent = child,
        px = cx,
        py = cy,
        center_angle = angles[[child]],
        width = max(42, width * 0.65),
        depth = depth + 1L
      )
    }
    invisible(NULL)
  }

  place_children(current, x, y, center_angle = 0, width = angle_width)

  nodes <- merge(tree$nodes, node_pos, by.x = "node_id", by.y = "node_id", all.x = TRUE, sort = FALSE)
  terminals <- ct_terminal_nodes(edges, tree$retained)
  nodes$is_terminal <- nodes$node_id %in% terminals
  nodes$sample_label <- ""
  for (terminal in terminals) {
    label <- paste(tree$present_samples[[terminal]], collapse = ", ")
    nodes$sample_label[nodes$node_id == terminal] <- label
  }

  list(nodes = nodes, edges = edge_pos)
}

#' @export
print.clonetree <- function(x, ...) {
  cat("CloneTree phylogeny\n")
  if (!is.null(x$input$subject)) {
    cat("Subject: ", x$input$subject, "\n", sep = "")
  }
  cat("Solution: ", x$solution, " of ", length(x$qc$solutions), "\n", sep = "")
  cat("Retained clusters: ", paste(x$retained, collapse = ", "), "\n", sep = "")
  if (length(x$excluded) > 0) {
    cat("Excluded clusters: ", paste(x$excluded, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
