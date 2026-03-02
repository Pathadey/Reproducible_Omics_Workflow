# ================================
# OMICS WORKFLOW v0.1 (STABLE)
# FUNCTIONS ONLY
# ================================

library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)

# ================================
# QIIME file parser
# ================================
process_qiime_file <- function(file_path) {
  
  lines <- readLines(file_path)
  
  if (length(lines) >= 2 && grepl("^#q2:types", lines[2])) {
    lines <- lines[-2]
  }
  
  tmp <- tempfile(fileext = ".tsv")
  writeLines(lines, tmp)
  
  df <- read_tsv(
    tmp,
    col_types = cols(.default = "c"),
    show_col_types = FALSE
  )
  
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  if (!"Taxon" %in% colnames(df)) {
    stop("Taxon column not found in ", file_path)
  }
  
  abundance_col <- tail(colnames(df), 1)
  
  out <- df %>%
    transmute(
      Taxon     = as.character(Taxon),
      Abundance = as.numeric(.data[[abundance_col]]),
      Sample    = sample_name
    )
  
  if (any(is.na(out$Abundance))) {
    stop("Non-numeric abundance values in ", file_path)
  }
  
  out
}

# ================================
# Constructor
# ================================
omics_init <- function(source = NULL) {
  
  omics <- list(
    raw_data = source,
    design   = NULL,
    otu      = NULL,
    alpha    = NULL,
    beta     = NULL,
    stats    = list(),
    meta     = list(
      created = Sys.time(),
      version = "v0.1"
    )
  )
  
  class(omics) <- "omics"
  omics
}

# ================================
# Import QIIME tables
# ================================
omics_from_qiime_long <- function(path, pattern = "\\.tsv$") {
  
  stopifnot(dir.exists(path))
  
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No QIIME .tsv files found")
  
  long <- map_dfr(files, process_qiime_file)
  
  otu <- xtabs(Abundance ~ Sample + Taxon, data = long)
  
  omics <- omics_init(source = long)
  omics$otu <- as.matrix(otu)
  
  omics
}

# ================================
# Metadata
# ================================
omics_make_metadata <- function(omics) {
  
  stopifnot(inherits(omics, "omics"))
  stopifnot(!is.null(omics$otu))
  
  data.frame(
    Sample = rownames(omics$otu),
    stringsAsFactors = FALSE
  )
}

omics_parse_sample_names <- function(
    metadata,
    sample_col = "Sample",
    sep = "_",
    into = c("Treatment", "Replicate")
) {
  
  parts <- strsplit(metadata[[sample_col]], sep)
  parts <- do.call(rbind, lapply(parts, `length<-`, length(into)))
  colnames(parts) <- into
  
  cbind(metadata, as.data.frame(parts, stringsAsFactors = FALSE))
}

# ================================
# Design
# ================================
omics_design <- function(
    omics,
    metadata,
    sample_col = "Sample",
    treatment_col = "Treatment"
) {
  
  stopifnot(inherits(omics, "omics"))
  stopifnot(!is.null(omics$otu))
  
  samples <- rownames(omics$otu)
  meta <- metadata[match(samples, metadata[[sample_col]]), ]
  
  if (length(unique(meta[[treatment_col]])) == nrow(meta)) {
    warning("No biological replicates detected; inferential statistics are descriptive.")
  }
  
  omics$design <- list(
    sample_col    = sample_col,
    treatment_col = treatment_col,
    metadata      = meta
  )
  
  omics
}

# ================================
# Alpha diversity
# ================================
omics_alpha <- function(omics) {
  
  otu <- omics$otu
  
  alpha <- data.frame(
    Sample   = rownames(otu),
    Shannon  = diversity(otu, "shannon"),
    Simpson  = diversity(otu, "simpson"),
    Observed = rowSums(otu > 0)
  )
  
  if (!is.null(omics$design)) {
    alpha <- merge(
      omics$design$metadata,
      alpha,
      by.x = omics$design$sample_col,
      by.y = "Sample"
    )
  }
  
  omics$alpha <- alpha
  omics
}

omics_alpha_stats <- function(
    omics,
    group_col = "Treatment"
) {
  
  stopifnot(inherits(omics, "omics"))
  stopifnot(!is.null(omics$alpha))
  stopifnot(!is.null(omics$design))
  
  alpha <- omics$alpha
  
  # Identify alpha metrics automatically
  non_metrics <- c(
    omics$design$sample_col,
    omics$design$treatment_col,
    "Replicate"
  )
  
  metrics <- setdiff(colnames(alpha), non_metrics)
  
  stats <- lapply(metrics, function(metric) {
    
    res <- kruskal.test(
      as.formula(paste(metric, "~", group_col)),
      data = alpha
    )
    
    data.frame(
      Metric    = metric,
      Test      = "Kruskal-Wallis",
      Statistic = unname(res$statistic),
      DF        = unname(res$parameter),
      P_value   = res$p.value
    )
  })
  
  omics$stats$alpha <- do.call(rbind, stats)
  omics
}

omics_alpha_summary <- function(
    omics,
    group_col = "Treatment"
) {
  
  stopifnot(!is.null(omics$alpha))
  
  alpha <- omics$alpha
  
  metrics <- c("Shannon", "Simpson", "Observed")
  
  summary <- aggregate(
    alpha[, metrics],
    by = list(Treatment = alpha[[group_col]]),
    FUN = function(x) c(
      mean = mean(x),
      sem  = sd(x) / sqrt(length(x))
    )
  )
  
  # Flatten columns
  out <- data.frame(Treatment = summary$Treatment)
  
  for (m in metrics) {
    out[[paste0(m, "_mean")]] <- summary[[m]][, "mean"]
    out[[paste0(m, "_SEM")]]  <- summary[[m]][, "sem"]
  }
  
  out
}

# ================================
# Beta diversity (STABLE version)
# ================================
omics_beta <- function(
    omics,
    distance = "bray",
    ordination = "NMDS"
) {
  
  otu  <- omics$otu
  meta <- omics$design$metadata
  trt  <- omics$design$treatment_col
  
  dist <- vegdist(otu, distance)
  
  scores <- NULL
  axis_var <- NULL
  stress <- NULL
  used <- ordination
  
  if (ordination == "NMDS") {
    ord <- suppressWarnings(metaMDS(dist, k = 2, trymax = 20))
    if (!is.na(ord$stress) && ord$stress > 1e-3 && nrow(otu) >= 8) {
      scores <- as.data.frame(ord$points)
      stress <- ord$stress
    } else {
      used <- "PCoA"
    }
  }
  
  if (used == "PCoA") {
    ord <- cmdscale(dist, k = 2, eig = TRUE)
    scores <- as.data.frame(ord$points)
    axis_var <- ord$eig[1:2]
  }
  
  colnames(scores) <- c("Axis1", "Axis2")
  scores$Sample <- rownames(scores)
  
  scores <- merge(
    meta,
    scores,
    by.x = omics$design$sample_col,
    by.y = "Sample"
  )
  
  meta$GroupVar <- meta[[trt]]
  permanova <- adonis2(dist ~ GroupVar, data = meta)
  
  omics$beta <- list(
    scores      = scores,
    permanova  = permanova,
    axis_var   = axis_var,
    stress     = stress,
    ordination = used
  )
  
  omics
}

# ================================
# Summary
# ================================
summary.omics <- function(object, ...) {
  
  cat("OMICS WORKFLOW SUMMARY\n")
  cat("Samples:", nrow(object$otu), "\n")
  cat("Features:", ncol(object$otu), "\n")
  
  if (!is.null(object$stats$alpha)) {
    print(object$stats$alpha)
  }
  
  if (!is.null(object$beta$permanova)) {
    print(object$beta$permanova)
  }
  
  invisible(object)
}
#Omics plot with centroids per replicate
omics_plot <- function(omics, group_col = "Treatment") {
  
  scores <- omics$beta$scores
  
  ggplot(
    scores,
    aes(x = Axis1, y = Axis2, colour = .data[[group_col]])
  ) +
    geom_point(size = 3) +
    theme_classic() +
    labs(
      title = paste("Ordination (", omics$beta$ordination, ")", sep = "")
    )
}

omics_abundance_table <- function(
    omics,
    rank = "Phylum",
    treatment_col = "Treatment",
    min_rel_abundance = 0.01
) {
  
  stopifnot(inherits(omics, "omics"))
  stopifnot(!is.null(omics$design))
  stopifnot(!is.null(omics$otu))
  
  # Long abundance table
  df <- as.data.frame(omics$otu)
  colnames(df) <- c("Sample", "Taxon", "Abundance")
  
  meta <- omics$design$metadata
  sample_col <- omics$design$sample_col
  
  df <- merge(meta, df, by.x = sample_col, by.y = "Sample")
  
  # Parse taxonomy
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  idx <- match(rank, ranks)
  
  if (is.na(idx)) {
    stop("Invalid rank. Choose from: ", paste(ranks, collapse = ", "))
  }
  
  tax_split <- strsplit(as.character(df$Taxon), ";")
  df$TaxRank <- sapply(tax_split, function(x) {
    if (length(x) >= idx) {
      gsub("^[a-z]__", "", trimws(x[idx]))
    } else {
      "Unclassified"
    }
  })
  
  # Aggregate by treatment
  agg <- aggregate(
    Abundance ~ .,
    data = df[, c(treatment_col, "TaxRank", "Abundance")],
    FUN = sum
  )
  
  colnames(agg)[1] <- "Treatment"
  
  # Relative abundance
  agg$Abundance <- ave(
    agg$Abundance,
    agg$Treatment,
    FUN = function(x) x / sum(x)
  )
  
  # Collapse rare taxa
  agg$TaxRank[agg$Abundance < min_rel_abundance] <- "Other"
  
  agg <- aggregate(
    Abundance ~ Treatment + TaxRank,
    data = agg,
    FUN = sum
  )
  
  # Wide table for Excel
  out <- reshape(
    agg,
    idvar = "Treatment",
    timevar = "TaxRank",
    direction = "wide"
  )
  
  colnames(out) <- gsub("^Abundance\\.", "", colnames(out))
  rownames(out) <- out$Treatment
  out$Treatment <- NULL
  
  out
}
plot_ordination <- function(
    omics,
    group_col = "Treatment"
) {
  
  scores <- omics$beta$scores
  
  ggplot(
    scores,
    aes(x = Axis1, y = Axis2, colour = .data[[group_col]])
  ) +
    geom_point(size = 3) +
    theme_classic() +
    labs(
      x = "Axis 1",
      y = "Axis 2",
      title = paste("Ordination (", omics$beta$ordination, ")", sep = "")
    )
}

