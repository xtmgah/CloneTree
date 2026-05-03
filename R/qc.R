#' Convert clustering results to CloneTree's internal format
#'
#' @description
#' Standardizes a long table of cancer cell fraction (CCF) estimates into the
#' matrix and metadata representation used by all CloneTree functions.
#'
#' @param data A data frame with one row per mutation cluster and sample.
#' @param subject Optional subject/tumor identifier to select when `data`
#'   contains more than one tumor.
#' @param sample_col,cluster_col,ccf_col Column names for sample ID, cluster
#'   ID, and CCF.
#' @param mutation_col Column with the number of substitutions/mutations in a
#'   cluster. If absent, each cluster is assigned length 1.
#' @param subject_col Optional subject ID column.
#' @param clone_col Optional column marking trunk clusters.
#' @param trunk_value Value in `clone_col` that marks a trunk cluster.
#' @param cluster_prefix Prefix added to purely numeric cluster IDs. Use `NULL`
#'   to keep cluster IDs unchanged.
#' @param ccf_scale One of `"auto"`, `"fraction"`, or `"percent"`.
#' @param snap_ccf Logical; if `TRUE`, values below `ccf_floor` become 0 and
#'   values above `ccf_ceiling` become 1.
#' @param ccf_floor,ccf_ceiling Thresholds used when `snap_ccf = TRUE`.
#' @param trunk_threshold When no `clone_col` is available, clusters with CCF
#'   at least this value in all samples are marked as trunk.
#' @param summarize_ccf How duplicate cluster/sample rows are summarized.
#'
#' @return A `clonetree_input` object.
#' @export
as_clonetree_input <- function(data,
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
  ct_standardize_input(
    data = data,
    subject = subject,
    sample_col = sample_col,
    cluster_col = cluster_col,
    ccf_col = ccf_col,
    mutation_col = mutation_col,
    subject_col = subject_col,
    clone_col = clone_col,
    trunk_value = trunk_value,
    cluster_prefix = cluster_prefix,
    ccf_scale = ccf_scale,
    snap_ccf = snap_ccf,
    ccf_floor = ccf_floor,
    ccf_ceiling = ccf_ceiling,
    trunk_threshold = trunk_threshold,
    summarize_ccf = summarize_ccf
  )
}

ct_pairwise_rules <- function(input,
                              clusters = rownames(input$ccf),
                              ccf_tolerance = 0.03,
                              crossing_tolerance = ccf_tolerance) {
  ccf <- input$ccf[clusters, , drop = FALSE]
  if (nrow(ccf) < 2) {
    return(data.frame(
      cluster_a = character(),
      cluster_b = character(),
      pigeonhole = logical(),
      crossing = logical(),
      incompatible = logical(),
      a_can_parent_b = logical(),
      b_can_parent_a = logical(),
      relation = character(),
      max_sum = numeric(),
      max_sum_sample = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- vector("list", choose(nrow(ccf), 2))
  k <- 1L
  sample_labels <- input$samples$sample_label[match(colnames(ccf), input$samples$sample_id)]
  for (i in seq_len(nrow(ccf) - 1L)) {
    for (j in (i + 1L):nrow(ccf)) {
      a <- ccf[i, ]
      b <- ccf[j, ]
      diff <- a - b
      sums <- a + b
      max_idx <- which.max(sums)
      crossing <- any(diff > crossing_tolerance) && any(diff < -crossing_tolerance)
      pigeonhole <- any(sums > 1 + ccf_tolerance)
      a_can_parent_b <- all(a + ccf_tolerance >= b) && any(a > b + ccf_tolerance)
      b_can_parent_a <- all(b + ccf_tolerance >= a) && any(b > a + ccf_tolerance)
      incompatible <- crossing && pigeonhole
      relation <- if (incompatible) {
        "incompatible"
      } else if (pigeonhole) {
        "nested_required"
      } else if (crossing) {
        "branch_required"
      } else if (a_can_parent_b || b_can_parent_a) {
        "nested_allowed"
      } else {
        "unconstrained"
      }
      out[[k]] <- data.frame(
        cluster_a = rownames(ccf)[i],
        cluster_b = rownames(ccf)[j],
        pigeonhole = pigeonhole,
        crossing = crossing,
        incompatible = incompatible,
        a_can_parent_b = a_can_parent_b,
        b_can_parent_a = b_can_parent_a,
        relation = relation,
        max_sum = unname(sums[max_idx]),
        max_sum_sample = sample_labels[max_idx],
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

ct_conflict_matrix <- function(nodes, pairwise) {
  mat <- matrix(FALSE, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  if (nrow(pairwise) > 0) {
    bad <- pairwise[pairwise$incompatible, , drop = FALSE]
    for (i in seq_len(nrow(bad))) {
      mat[bad$cluster_a[i], bad$cluster_b[i]] <- TRUE
      mat[bad$cluster_b[i], bad$cluster_a[i]] <- TRUE
    }
  }
  mat
}

ct_best_independent_sets <- function(nodes,
                                     conflicts,
                                     weights,
                                     locked = character(),
                                     max_solutions = 1000) {
  nodes <- unique(as.character(nodes))
  locked <- intersect(unique(as.character(locked)), nodes)
  weights <- weights[nodes]
  weights[is.na(weights)] <- 1
  if (length(locked) > 1) {
    locked_conflicts <- conflicts[locked, locked, drop = FALSE]
    diag(locked_conflicts) <- FALSE
    if (any(locked_conflicts)) {
      stop("Locked clusters are mutually incompatible: ",
           paste(locked, collapse = ", "),
           call. = FALSE)
    }
  }

  locked_neighbors <- character()
  if (length(locked) > 0) {
    locked_neighbors <- names(which(colSums(conflicts[locked, , drop = FALSE]) > 0))
  }
  optional <- setdiff(nodes, union(locked, locked_neighbors))
  degree <- rowSums(conflicts[optional, optional, drop = FALSE])
  optional <- optional[order(-degree, -weights[optional], optional)]

  best_score <- -Inf
  best <- list()
  truncated <- FALSE

  add_solution <- function(sol, score) {
    sol <- nodes[nodes %in% sol]
    if (score > best_score + .Machine$double.eps) {
      best_score <<- score
      best <<- list(sol)
      truncated <<- FALSE
    } else if (abs(score - best_score) <= .Machine$double.eps) {
      if (length(best) < max_solutions) {
        best[[length(best) + 1L]] <<- sol
      } else {
        truncated <<- TRUE
      }
    }
  }

  search <- function(included, candidates, score) {
    if (length(candidates) == 0) {
      add_solution(included, score)
      return(invisible(NULL))
    }
    upper <- score + sum(weights[candidates])
    if (upper < best_score - .Machine$double.eps) {
      return(invisible(NULL))
    }

    v <- candidates[[1]]
    remaining <- candidates[-1]
    neighbors <- names(which(conflicts[v, ]))
    search(
      included = c(included, v),
      candidates = setdiff(remaining, neighbors),
      score = score + weights[[v]]
    )
    search(
      included = included,
      candidates = remaining,
      score = score
    )
    invisible(NULL)
  }

  search(
    included = locked,
    candidates = optional,
    score = sum(weights[locked])
  )

  if (length(best) == 0) {
    best <- list(locked)
    best_score <- sum(weights[locked])
  }

  list(
    solutions = best,
    score = best_score,
    truncated = truncated,
    forced_excluded = locked_neighbors
  )
}

#' Check clone clusters against pigeonhole and crossing rules
#'
#' @description
#' `clone_qc()` evaluates pairwise CCF relationships between mutation clusters.
#' A pair is considered incompatible when the two clusters cross across samples
#' and also violate the pigeonhole rule in at least one sample. In that case the
#' two clusters cannot both be kept under the usual evolutionary assumptions, so
#' CloneTree enumerates all best-scoring compatible retained-cluster sets.
#'
#' @param data A clustering result table or a `clonetree_input` object.
#' @param subject Optional tumor/subject ID.
#' @param remove_clusters Optional cluster IDs to remove before QC. IDs can be
#'   either original IDs (for example `"11"`) or prefixed IDs (for example
#'   `"C11"`).
#' @param presence_threshold Clusters with maximum CCF below this value are
#'   treated as absent.
#' @param ccf_tolerance Tolerance for the pigeonhole rule. Two independent
#'   clusters are incompatible with separate branches if their summed CCF is
#'   greater than `1 + ccf_tolerance` in any sample.
#' @param crossing_tolerance Tolerance for the crossing rule. Defaults to
#'   `ccf_tolerance`.
#' @param objective How to choose among possible exclusions. `"max_mutations"`
#'   keeps compatible solutions with the largest retained mutation count;
#'   `"max_clusters"` keeps the largest number of clusters.
#' @param max_solutions Maximum number of equally best compatible solutions to
#'   keep.
#' @param lock_trunk If `TRUE`, trunk clusters are never excluded by the solver.
#' @param locked_clusters Additional clusters that must be retained.
#' @param ... Passed to [as_clonetree_input()] when `data` is a raw data frame.
#'
#' @return A `clonetree_qc` object containing the standardized input, pairwise
#'   rule table, and all retained-cluster solutions.
#' @export
clone_qc <- function(data,
                     subject = NULL,
                     remove_clusters = NULL,
                     presence_threshold = 0.02,
                     ccf_tolerance = 0.03,
                     crossing_tolerance = ccf_tolerance,
                     objective = c("max_mutations", "max_clusters"),
                     max_solutions = 1000,
                     lock_trunk = TRUE,
                     locked_clusters = NULL,
                     ...) {
  objective <- match.arg(objective)
  input <- ct_get_input(data, subject = subject, ...)
  meta <- input$meta
  active <- meta$cluster_id[meta$max_ccf >= presence_threshold]
  manual_removed <- ct_normalize_requested_clusters(remove_clusters, meta, input$cluster_prefix)
  active <- setdiff(active, manual_removed)
  if (length(active) == 0) {
    stop("No active clusters remain after filtering/removal.", call. = FALSE)
  }

  pairwise <- ct_pairwise_rules(
    input = input,
    clusters = active,
    ccf_tolerance = ccf_tolerance,
    crossing_tolerance = crossing_tolerance
  )
  conflicts <- ct_conflict_matrix(active, pairwise)
  weights <- if (objective == "max_clusters") {
    stats::setNames(rep(1, nrow(meta)), meta$cluster_id)
  } else {
    stats::setNames(meta$mutations, meta$cluster_id)
  }

  locked <- ct_normalize_requested_clusters(locked_clusters, meta, input$cluster_prefix)
  if (lock_trunk) {
    locked <- union(locked, meta$cluster_id[meta$is_trunk])
  }
  locked <- intersect(locked, active)

  best <- ct_best_independent_sets(
    nodes = active,
    conflicts = conflicts,
    weights = weights,
    locked = locked,
    max_solutions = max_solutions
  )

  solutions <- lapply(best$solutions, function(sol) {
    sol <- active[active %in% sol]
    list(
      retained = sol,
      excluded = setdiff(active, sol),
      score = sum(weights[sol])
    )
  })

  out <- list(
    input = input,
    active_clusters = active,
    manual_removed = manual_removed,
    pairwise = pairwise,
    conflicts = conflicts,
    solutions = solutions,
    objective = objective,
    ccf_tolerance = ccf_tolerance,
    crossing_tolerance = crossing_tolerance,
    presence_threshold = presence_threshold,
    max_solutions = max_solutions,
    truncated = best$truncated,
    locked_clusters = locked
  )
  class(out) <- "clonetree_qc"
  out
}

#' Extract retained-cluster solutions from a QC result
#'
#' @param qc A `clonetree_qc` object.
#' @return A data frame with one row per solution and comma-separated retained
#'   and excluded cluster IDs.
#' @export
qc_solutions <- function(qc) {
  if (!inherits(qc, "clonetree_qc")) {
    stop("`qc` must be a clonetree_qc object.", call. = FALSE)
  }
  data.frame(
    solution = seq_along(qc$solutions),
    score = vapply(qc$solutions, `[[`, numeric(1), "score"),
    retained = vapply(qc$solutions, function(x) paste(x$retained, collapse = ", "), character(1)),
    excluded = vapply(qc$solutions, function(x) paste(x$excluded, collapse = ", "), character(1)),
    stringsAsFactors = FALSE
  )
}

#' @export
print.clonetree_qc <- function(x, ...) {
  n_conflicts <- sum(x$pairwise$incompatible)
  cat("CloneTree QC\n")
  if (!is.null(x$input$subject)) {
    cat("Subject: ", x$input$subject, "\n", sep = "")
  }
  cat("Active clusters: ", length(x$active_clusters), "\n", sep = "")
  cat("Incompatible pairs: ", n_conflicts, "\n", sep = "")
  cat("Solutions retained: ", length(x$solutions), "\n", sep = "")
  if (x$truncated) {
    cat("Note: solution list was truncated at max_solutions.\n")
  }
  invisible(x)
}
