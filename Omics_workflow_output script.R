# ================================
# EXPORT ALL RESULTS
# ================================

stopifnot(exists("omics"))

# ----------------
# 1) ALPHA DIVERSITY
# ----------------
dir.create("results/alpha", recursive = TRUE, showWarnings = FALSE)

write.csv(
  omics$alpha,
  "results/alpha/alpha_diversity_values.csv",
  row.names = FALSE
)

write.csv(
  omics$stats$alpha,
  "results/alpha/alpha_diversity_stats.csv",
  row.names = FALSE
)

# Optional summary (mean ± SEM)
alpha_summary <- omics_alpha_summary(omics)

write.csv(
  alpha_summary,
  "results/alpha/alpha_diversity_mean_sem.csv",
  row.names = FALSE
)

# ----------------
# 2) ORDINATION (PCoA / NMDS)
# ----------------
dir.create("results/ordination", recursive = TRUE, showWarnings = FALSE)

# Ordination coordinates (THIS IS WHAT YOU PLOT)
write.csv(
  omics$beta$scores,
  "results/ordination/ordination_scores.csv",
  row.names = FALSE
)

# PERMANOVA
write.csv(
  as.data.frame(omics$beta$permanova),
  "results/ordination/ordination_permanova.csv"
)

# PCoA axis variance (if applicable)
if (!is.null(omics$beta$axis_var)) {
  write.csv(
    data.frame(
      Axis = c("Axis1", "Axis2"),
      Variance = omics$beta$axis_var
    ),
    "results/ordination/ordination_axis_variance.csv",
    row.names = FALSE
  )
}

# NMDS stress (if applicable)
if (!is.null(omics$beta$stress)) {
  write.csv(
    data.frame(Stress = omics$beta$stress),
    "results/ordination/ordination_nmds_stress.csv",
    row.names = FALSE
  )
}

# ----------------
# 3) TAXONOMIC ABUNDANCE TABLES
# ----------------
dir.create("results/abundance_tables", recursive = TRUE, showWarnings = FALSE)

taxonomic_ranks <- c(
  "Kingdom","Phylum","Class","Order","Family","Genus","Species"
)

for (rk in taxonomic_ranks) {
  
  message("Exporting abundance table at ", rk, " level")
  
  tab <- try(
    omics_abundance_table(
      omics,
      rank = rk,
      min_rel_abundance = 0.02
    ),
    silent = TRUE
  )
  
  if (inherits(tab, "try-error")) {
    warning("Skipping ", rk, " (no data)")
    next
  }
  
  write.csv(
    tab,
    file = paste0(
      "results/abundance_tables/abundance_",
      tolower(rk),
      ".csv"
    )
  )
}
