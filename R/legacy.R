#' Legacy oval CCF QC plot wrapper
#'
#' Compatibility wrapper for the original draft function name. New code should
#' prefer [clone_qc()] followed by [plot_clone_ccf()].
#'
#' @param cluster_data_ind Original cluster table.
#' @param subj Subject ID.
#' @param remove_info Cluster IDs to remove before QC.
#' @param clone_fill_colors Optional colors.
#' @param oval_width,oval_heigh Ellipse radii multipliers from the original
#'   draft function.
#' @param filename Optional output filename.
#' @param fwidth,fheight Optional output dimensions.
#' @param text_size Base text size for the oval plot.
#' @param ... Ignored legacy plotting arguments.
#'
#' @return A `ggplot` object.
#' @export
oval_clone_QC <- function(cluster_data_ind,
                          subj,
                          remove_info = NULL,
                          clone_fill_colors = NULL,
                          oval_width = 6,
                          oval_heigh = 10,
                          filename = NULL,
                          fwidth = NULL,
                          fheight = NULL,
                          text_size = 11,
                          ...) {
  qc <- clone_qc(cluster_data_ind, subject = subj, remove_clusters = remove_info)
  plot_clone_ccf(
    qc,
    solution = 1,
    show_values = TRUE,
    clonesum = FALSE,
    include_excluded = TRUE,
    solution_position = "subtitle",
    palette = clone_fill_colors %||% clonetree_palette(),
    text_size = text_size,
    oval_width = oval_width,
    oval_height = oval_heigh,
    filename = filename,
    width = fwidth,
    height = fheight
  )
}

#' Legacy oval CCF plot wrapper
#'
#' Compatibility wrapper for the original draft function name. New code should
#' prefer [plot_clone_ccf()].
#'
#' @param cluster_data_ind Original cluster table.
#' @param subj Subject ID.
#' @param remove_info Cluster IDs to remove before QC.
#' @param clone_fill_colors Optional colors.
#' @param oval_width,oval_heigh Ellipse radii multipliers from the original
#'   draft function.
#' @param clonesum If `TRUE`, add a clone-summary row.
#' @param filename Optional output filename.
#' @param fwidth,fheight Optional output dimensions.
#' @param text_size Base text size for the oval plot.
#' @param ... Ignored legacy plotting arguments.
#'
#' @return A `ggplot` object or saved filename invisibly.
#' @export
oval_clone_plot <- function(cluster_data_ind,
                            subj,
                            remove_info = NULL,
                            clone_fill_colors = NULL,
                            oval_width = 6,
                            oval_heigh = 10,
                            clonesum = TRUE,
                            filename = NULL,
                            fwidth = NULL,
                            fheight = NULL,
                            text_size = 11,
                            ...) {
  qc <- clone_qc(cluster_data_ind, subject = subj, remove_clusters = remove_info)
  plot_clone_ccf(
    qc,
    solution = 1,
    show_values = FALSE,
    clonesum = clonesum,
    include_excluded = FALSE,
    solution_position = "caption",
    subclone_outline = TRUE,
    summary_subclone_alpha = 0.75,
    add_group_spacing = TRUE,
    palette = clone_fill_colors %||% clonetree_palette(),
    text_size = text_size,
    oval_width = oval_width,
    oval_height = oval_heigh,
    filename = filename,
    width = fwidth,
    height = fheight
  )
}

#' Legacy tree plot wrapper
#'
#' Compatibility wrapper for the original draft function name. New code should
#' prefer [clone_qc()], [infer_clone_tree()], and [plot_clone_tree()].
#'
#' @param cluster_data_ind Original cluster table.
#' @param subj Subject ID.
#' @param cluster_data_info Optional event annotation table.
#' @param remove_info Cluster IDs to remove before QC.
#' @param clone_fill_colors Optional colors.
#' @param filename Optional output filename.
#' @param fwidth,fheight Optional output dimensions.
#' @param ... Ignored legacy plotting arguments.
#'
#' @return A `ggplot` object or saved filename invisibly.
#' @export
tree_base_plot <- function(cluster_data_ind,
                           subj,
                           cluster_data_info = NULL,
                           remove_info = NULL,
                           clone_fill_colors = NULL,
                           filename = NULL,
                           fwidth = NULL,
                           fheight = NULL,
                           ...) {
  qc <- clone_qc(cluster_data_ind, subject = subj, remove_clusters = remove_info)
  tree <- infer_clone_tree(qc, solution = 1)
  plot_clone_tree(
    tree,
    events = cluster_data_info,
    palette = clone_fill_colors %||% clonetree_palette(),
    filename = filename,
    width = fwidth,
    height = fheight
  )
}
