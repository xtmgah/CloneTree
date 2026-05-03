#' Example multi-region mutation cluster CCF table
#'
#' A small example dataset from the original CloneTree draft. Each row is a
#' mutation cluster in one tumor region/sample.
#'
#' @format A data frame with 40 rows and 7 columns:
#' \describe{
#'   \item{cluster.no}{Mutation cluster identifier.}
#'   \item{Clone}{`"Y"` for trunk clone clusters and `"N"` for subclonal clusters.}
#'   \item{no.of.mutations}{Number of mutations assigned to the cluster.}
#'   \item{meanccf}{Mean CCF across sampled tumor regions.}
#'   \item{sampleID}{Sample/region identifier.}
#'   \item{CCF}{Cancer cell fraction, represented as a fraction from 0 to 1.}
#'   \item{subject}{Subject/tumor identifier.}
#' }
#'
#' @source `CloneTree_Original_Package/cluster_data_ind2.RData`
#' @docType data
#' @usage data(cluster_data_ind2)
"cluster_data_ind2"
