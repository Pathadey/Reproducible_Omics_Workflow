omics_compute_centroids <- function(
    omics,
    group_col = "Treatment"
) {
  
  stopifnot(!is.null(omics$beta))
  stopifnot(!is.null(omics$beta$scores))
  stopifnot(group_col %in% colnames(omics$beta$scores))
  
  scores <- omics$beta$scores
  
  centroids <- aggregate(
    cbind(Axis1, Axis2) ~ scores[[group_col]],
    data = scores,
    FUN = mean
  )
  
  colnames(centroids)[1] <- group_col
  centroids
}


omics_plot_ordination <- function(
    omics,
    group_col = "Treatment",
    show_centroids = TRUE
) {
  
  stopifnot(!is.null(omics$beta))
  stopifnot(!is.null(omics$beta$scores))
  
  scores <- omics$beta$scores
  
  p <- ggplot(
    scores,
    aes(x = Axis1, y = Axis2, colour = .data[[group_col]])
  ) +
    geom_point(size = 3, alpha = 0.8) +
    theme_classic() +
    labs(
      title = paste("Ordination (", omics$beta$ordination, ")", sep = ""),
      x = "Axis 1",
      y = "Axis 2"
    )
  
  if (show_centroids) {
    
    if (is.null(omics$beta$centroids)) {
      omics$beta$centroids <- omics_compute_centroids(omics, group_col)
    }
    
    p <- p +
      geom_point(
        data = omics$beta$centroids,
        aes(x = Axis1, y = Axis2),
        shape = 4,
        size = 5,
        stroke = 1.2,
        colour = "black"
      )
  }
  
  p
}
