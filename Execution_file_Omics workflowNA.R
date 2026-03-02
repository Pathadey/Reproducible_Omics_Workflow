# ================================
# RUN OMICS ANALYSIS
# ================================

source("Omics_WorkflowNA.R")

# ---- Import data ----
omics <- omics_from_qiime_long("Data")

# ---- Metadata ----
meta <- omics_make_metadata(omics)
meta <- omics_parse_sample_names(meta, into = c("Treatment", "Replicate"))

omics <- omics_design(omics, meta)

# ---- Alpha diversity ----
omics <- omics_alpha(omics)
omics <- omics_alpha_stats(omics)

# ---- Beta diversity ----
omics <- omics_beta(omics)

summary(omics)
plot_ordination(omics) #per replicate
source("Omic_helper.R")
omics$beta$centroids <- omics_compute_centroids(
  omics,
  group_col = "Treatment"
)
omics$beta$centroids
p <- omics_plot_ordination(omics)
p

scores <- omics$beta$scores
scores$Type <- "Sample"

centroids <- omics$beta$centroids

p <- ggplot(
  centroids,
  aes(
    x = Axis1,
    y = Axis2,
    colour = Treatment
  )
) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    title = paste("Ordination (", omics$beta$ordination, ")", sep = ""),
    x = "Axis 1",
    y = "Axis 2"
  )

p

cat("Exporting results...\n")

source("Omics_workflow_output script.R")

cat("All results saved in /results folder.\n")

# Make sure results folder exists
dir.create("results", showWarnings = FALSE)

# Save plot
ggsave("results/nmds_plot.png", plot = p, width = 8, height = 6, dpi = 300)


