library(CloneTree)

data(cluster_data_ind2)

qc <- clone_qc(cluster_data_ind2, subject = "IGC-11-1130")
print(qc_solutions(qc))

plot_clone_qc(qc)
plot_qc_rules_guide()
plot_clone_ccf(qc, solution = 1)
oval_clone_QC(cluster_data_ind2, "IGC-11-1130")
oval_clone_plot(cluster_data_ind2, "IGC-11-1130")

events <- data.frame(
  subject = "IGC-11-1130",
  cluster.no = c("C1", "C8", "C10"),
  info = c("Truncal driver\nchr7 ASCNA", "chr2 ASCNA", "MET p.X123")
)

tree <- infer_clone_tree(qc, solution = 1)
plot_clone_tree(tree, events = events)
