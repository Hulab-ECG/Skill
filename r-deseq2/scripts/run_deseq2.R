#!/usr/bin/env Rscript
# Standardized DESeq2 wrapper for differential expression analysis
#
# Required inputs:
#   --count-matrix  gene-level raw count matrix CSV (gene_id in first column)
#   --meta          sample metadata TSV with sample, group columns
#   --formula       condition formula in DESeq2 syntax, e.g. ~group
#   --contrasts     one or more contrasts as "numerator,denominator[,label]"
#
# Optional:
#   --outdir       output directory (default: deseq2_results)
#   --fdr          adjusted p-value cutoff (default: 0.05)
#   --log2fc       absolute log2 fold change cutoff (default: 1)
#   --min-count    low-count filter: minimum count value (default: 946)
#   --min-samples  low-count filter: minimum number of samples above min-count (default: 2)
#   --genes        candidate gene list file; if provided, UGT-like subsets are also written
#   --alpha        DESeq2 results alpha if different from --fdr
#   --shrink       apply apeglm/normal/shrinkage log2FC (NULL/none|normal|apeglm|ashr)
#   --lfc-shrink   alias for --shrink
#
# Example:
#   Rscript run_deseq2.R \
#     --count-matrix gene_count_matrix.csv \
#     --meta sample_meta.tsv \
#     --formula ~condition \
#     --contrasts blue,control,blue_vs_control \
#     --contrasts red,control,red_vs_control \
#     --contrasts red,blue,red_vs_blue \
#     --genes ugt_candidates.txt \
#     --outdir deseq2_results

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is not installed in the active R library.\n",
         "Install it with BiocManager::install(\"DESeq2\")")
  }
  library(DESeq2)
})

# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- list(
    count_matrix = NULL,
    meta = NULL,
    formula = "~group",
    contrasts = list(),
    outdir = "deseq2_results",
    fdr = 0.05,
    log2fc = 1.0,
    min_count = 10,
    min_samples = 2,
    genes = NULL,
    alpha = NULL,
    shrink = "none"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    take_val <- key %in% c(
      "--count-matrix", "--meta", "--formula", "--contrasts", "--outdir",
      "--fdr", "--log2fc", "--min-count", "--min-samples", "--genes",
      "--alpha", "--shrink", "--lfc-shrink"
    )
    if (!take_val) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    switch(key,
      "--count-matrix" = { opts$count_matrix <- val; i <- i + 2 },
      "--meta"         = { opts$meta <- val; i <- i + 2 },
      "--formula"      = { opts$formula <- val; i <- i + 2 },
      "--contrasts"    = {
        # Accept either one or several contrasts in the same argument, e.g.
        #   --contrasts blue,control blue,control,blue_vs_control
        #   --contrasts blue,control red,control
        toks <- strsplit(val, "[[:space:],]+")[[1]]
        toks <- trimws(toks)
        toks <- toks[nchar(toks) > 0]
        j <- 1
        while (j <= length(toks)) {
          if (j + 2 <= length(toks) &&
              grepl("_vs_", paste(toks[(j + 2):length(toks)], collapse = " ")[1])) {
            # Three-part triplet: numerator, denominator, label
            opts$contrasts[[length(opts$contrasts) + 1]] <- list(
              numerator = toks[j],
              denominator = toks[j + 1],
              label = toks[j + 2]
            )
            j <- j + 3
          } else if (j + 1 <= length(toks)) {
            opts$contrasts[[length(opts$contrasts) + 1]] <- list(
              numerator = toks[j],
              denominator = toks[j + 1],
              label = paste0(toks[j], "_vs_", toks[j + 1])
            )
            j <- j + 2
          } else {
            stop("Each --contrasts must be 'numerator,denominator' or ",
                 "'numerator,denominator,label'", call. = FALSE)
          }
        }
        i <- i + 2
      },
      "--outdir"      = { opts$outdir <- val; i <- i + 2 },
      "--fdr"         = { opts$fdr <- as.numeric(val); i <- i + 2 },
      "--log2fc"      = { opts$log2fc <- as.numeric(val); i <- i + 2 },
      "--min-count"   = { opts$min_count <- as.numeric(val); i <- i + 2 },
      "--min-samples" = { opts$min_samples <- as.integer(val); i <- i + 2 },
      "--genes"       = { opts$genes <- val; i <- i + 2 },
      "--alpha"       = { opts$alpha <- as.numeric(val); i <- i + 2 },
      "--shrink"      = { opts$shrink <- tolower(val); i <- i + 2 },
      "--lfc-shrink"  = { opts$shrink <- tolower(val); i <- i + 2 }
    )
  }

  if (is.null(opts$count_matrix) || !file.exists(opts$count_matrix)) {
    stop("Count matrix not found: ", opts$count_matrix, call. = FALSE)
  }
  if (is.null(opts$meta)) {
    stop("Missing required --meta", call. = FALSE)
  }
  if (!file.exists(opts$meta)) {
    stop("Metadata not found: ", opts$meta, call. = FALSE)
  }
  if (length(opts$contrasts) == 0) {
    stop("At least one --contrasts is required", call. = FALSE)
  }
  opts$alpha <- if (is.null(opts$alpha)) opts$fdr else opts$alpha
  opts
}

# ------------------------------------------------------------------
# Input helpers
# ------------------------------------------------------------------
read_counts <- function(path) {
  cm <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(cm) < 2) {
    stop("Count matrix has no sample columns: ", path, call. = FALSE)
  }
  gid_col <- names(cm)[1]
  counts <- as.matrix(cm[, -1, drop = FALSE])
  rownames(counts) <- cm[[gid_col]]
  mode(counts) <- "numeric"
  counts[is.na(counts)] <- 0
  counts <- round(counts)
  list(counts = counts, samples = setdiff(names(cm), gid_col))
}

read_metadata <- function(path, count_samples) {
  m <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% colnames(m)) {
    stop("Metadata must contain a 'sample' column.", call. = FALSE)
  }
  rownames(m) <- m$sample
  m <- m[count_samples, , drop = FALSE]
  m
}

read_gene_list <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(character(0))
  }
  ids <- trimws(readLines(path, warn = FALSE))
  ids <- ids[nchar(ids) > 0]
  unique(ids)
}

# ------------------------------------------------------------------
# Core analysis function
# ------------------------------------------------------------------
run_standard_deseq2 <- function(count_matrix_path,
                                meta_path,
                                formula = "~group",
                                contrasts = list(),
                                outdir = "deseq2_results",
                                fdr = 0.05,
                                log2fc = 1.0,
                                min_count = 10,
                                min_samples = 2,
                                gene_list_path = NULL,
                                alpha = NULL,
                                shrink = "none") {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  alpha <- if (is.null(alpha)) fdr else alpha

  # load inputs
  cm_info <- read_counts(count_matrix_path)
  counts <- cm_info$counts
  sample_names <- cm_info$samples
  meta <- read_metadata(meta_path, sample_names)
  gene_ids <- read_gene_list(gene_list_path)

  cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")
  cat("Samples used:", paste(sample_names, collapse = ", "), "\n")
  print(meta)

  # validate formula terms are present in metadata
  if (!"group" %in% colnames(meta) && formula == "~group") {
    stop("Metadata does not contain a 'group' column. ",
         "Use --formula to set the DESeq2 design.", call. = FALSE)
  }
  formula_obj <- as.formula(formula)
  all_vars <- all.vars(formula_obj)
  missing_vars <- setdiff(all_vars, colnames(meta))
  if (length(missing_vars) > 0) {
    stop("Metadata missing column(s) used in --formula: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  # make sure covariates/factors are factors, except first condition-like variable
  for (v in all_vars) {
    if (!is.factor(meta[[v]])) meta[[v]] <- factor(meta[[v]])
  }

  # The last design variable is treated as the main grouping factor when
  # contrasts are written as numerator,denominator[,label]. This matches the
  # common DESeq2 pattern of a single condition factor in the design, e.g.
  # --formula ~group or --formula ~batch+group.
  condition_var <- tail(all_vars, 1)

  col_data <- data.frame(meta, row.names = rownames(meta), check.names = FALSE)
  col_data <- col_data[sample_names, , drop = FALSE]

  # low-count filter
  keep <- rowSums(counts >= min_count) >= min_samples
  cat("Genes passing low-count filter", min_count, "in >=", min_samples,
      "samples:", sum(keep), "\n")
  counts <- counts[keep, , drop = FALSE]
  if (nrow(counts) == 0) {
    stop("No genes passed the low-count filter", call. = FALSE)
  }

  # DESeq
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data,
                                design = formula_obj)
  dds <- DESeq(dds, quiet = TRUE)

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts <- data.frame(gene_id = rownames(norm_counts), norm_counts,
                            check.names = FALSE)
  write.csv(norm_counts, file.path(outdir, "normalized_counts.csv"),
            row.names = FALSE)
  cat("Normalized counts written:", file.path(outdir, "normalized_counts.csv"), "\n")

  # candidate gene subset naming
  gene_label <- "cand"
  if (length(gene_ids) > 0) {
    gene_ids <- intersect(gene_ids, rownames(counts))
    cat("Candidate gene IDs present in count matrix:", length(gene_ids), "\n")
  }

  run_comparison <- function(contrast, shrink_method) {
    numerator <- contrast$numerator
    denominator <- contrast$denominator
    label <- contrast$label

    res <- results(dds, contrast = c(condition_var, numerator, denominator),
                   alpha = alpha)

    if (shrink_method %in% c("normal", "apeglm", "ashr")) {
      if (requireNamespace("apeglm", quietly = TRUE) || shrink_method != "apeglm") {
        tryCatch({
          res <- lfcShrink(dds,
                           contrast = c(condition_var, numerator, denominator),
                           res = res, type = shrink_method)
        }, error = function(e) {
          warning("lfcShrink failed (", conditionMessage(e), "); keeping unshrunk LFC.")
        })
      } else {
        warning("apeglm not installed; keeping unshrunk LFC.")
      }
    }

    df <- as.data.frame(res)
    df$gene_id <- rownames(df)
    df <- df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                 "stat", "pvalue", "padj")]
    df$significant <- !is.na(df$padj) & df$padj < fdr &
      abs(df$log2FoldChange) >= log2fc
    df <- df[order(df$padj, na.last = TRUE), ]

    write.csv(df, file.path(outdir, paste0("deseq2_", label, "_all.csv")),
              row.names = FALSE)

    if (length(gene_ids) > 0) {
      sub <- df[df$gene_id %in% gene_ids, , drop = FALSE]
      write.csv(sub, file.path(outdir, paste0("deseq2_", label, "_", gene_label, ".csv")),
                row.names = FALSE)
    }

    n_sig <- sum(df$significant, na.rm = TRUE)
    cat(label, ": ", n_sig, " significant genes (FDR<", fdr,
        ", |log2FC|>=", log2fc, ")\n", sep = "")
    list(df = df, label = label)
  }

  comparisons <- lapply(contrasts, run_comparison, shrink_method = shrink)
  names(comparisons) <- vapply(contrasts, function(x) x$label, "")

  # candidate gene union summary
  if (length(gene_ids) > 0) {
    summary_df <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    for (comp in comparisons) {
      df <- comp$df
      sub <- df[df$gene_id %in% gene_ids,
                c("gene_id", "log2FoldChange", "padj", "significant")]
      colnames(sub) <- c("gene_id",
                         paste0("log2FC_", comp$label),
                         paste0("padj_", comp$label),
                         paste0("sig_", comp$label))
      summary_df <- merge(summary_df, sub, by = "gene_id", all.x = TRUE)
    }
    sig_cols <- grep("^sig_", names(summary_df), value = TRUE)
    summary_df$n_significant <- rowSums(summary_df[, sig_cols, drop = FALSE],
                                        na.rm = TRUE)
    padj_cols <- grep("^padj_", names(summary_df), value = TRUE)
    summary_df$min_padj <- apply(summary_df[, padj_cols, drop = FALSE], 1, function(x) {
      v <- x[!is.na(x)]
      if (length(v) > 0) min(v) else NA_real_
    })
    summary_df <- summary_df[summary_df$n_significant > 0, , drop = FALSE]
    summary_df <- summary_df[order(-summary_df$n_significant,
                                   summary_df$min_padj), , drop = FALSE]
    write.csv(summary_df,
              file.path(outdir, paste0("significant_", gene_label, "_union.csv")),
              row.names = FALSE)
  }

  # text report
  rpt <- file.path(outdir, "de_summary.txt")
  sink(rpt)
  cat("DESeq2 differential expression summary\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Design:", formula, "\n")
  cat("Thresholds: FDR<", fdr, ", |log2FC|>=", log2fc,
      ", min count=", min_count, " in ", min_samples, " samples\n\n", sep = "")

  for (comp in comparisons) {
    df <- comp$df
    nm <- comp$label
    all_sig <- sum(df$significant, na.rm = TRUE)
    up <- sum(df$significant & df$log2FoldChange > 0, na.rm = TRUE)
    down <- sum(df$significant & df$log2FoldChange < 0, na.rm = TRUE)
    cat(nm, ": all DE genes =", all_sig, " (up=", up, ", down=", down, ")\n",
        sep = "")
    if (length(gene_ids) > 0) {
      ugt <- df[df$gene_id %in% gene_ids, , drop = FALSE]
      ugt_sig <- ugt[!is.na(ugt$padj) & ugt$padj < fdr &
                       abs(ugt$log2FoldChange) >= log2fc, , drop = FALSE]
      cat("  candidate DE =", nrow(ugt_sig),
          " (up=", sum(ugt_sig$log2FoldChange > 0),
          ", down=", sum(ugt_sig$log2FoldChange < 0), ")\n", sep = "")
      if (nrow(ugt_sig) > 0) {
        ugt_sig <- ugt_sig[order(ugt_sig$padj), , drop = FALSE]
        cat("  Top candidates:\n")
        for (i in seq_len(min(20, nrow(ugt_sig)))) {
          cat("    ", ugt_sig$gene_id[i], " log2FC=",
              round(ugt_sig$log2FoldChange[i], 3), " padj=",
              format(ugt_sig$padj[i], scientific = TRUE, digits = 3), "\n",
              sep = "")
        }
      }
      cat("\n")
    }
  }
  sink()
  cat("Summary report:", rpt, "\n")

  # plots for each contrast
  plot_one <- function(df, nm) {
    p <- ifelse(is.na(df$pvalue), 1, df$pvalue)
    is_sig <- !is.na(df$significant) & df$significant

    pdf(file.path(outdir, paste0("volcano_", nm, ".pdf")),
        width = 10, height = 8.9)
    plot(df$log2FoldChange, -log10(p), pch = 20,
         col = ifelse(is_sig, "firebrick", "grey70"),
         xlab = "log2 fold change", ylab = "-log10 p-value",
         main = paste("Volcano:", gsub("_", " vs ", nm)))
    abline(h = -log10(fdr), lty = 2, col = "grey30")
    abline(v = c(-log2fc, log2fc), lty = 2, col = "grey30")
    if (length(gene_ids) > 0) {
      ugt_idx <- which(df$gene_id %in% gene_ids)
      if (length(ugt_idx) > 0) {
        ord <- order(df$padj[ugt_idx], na.last = TRUE)
        top_ugt <- head(ugt_idx[ord], 10)
        text(df$log2FoldChange[top_ugt], -log10(p[top_ugt]),
             labels = df$gene_id[top_ugt], cex = 0.55, pos = 3)
      }
    }
    legend("topright", legend = c("significant", "candidates labelled"),
           pch = c(20, NA), col = c("firebrick", "black"),
           text.col = c("firebrick", "black"), cex = 0.7, bty = "n")
    dev.off()

    pdf(file.path(outdir, paste0("ma_", nm, ".pdf")),
        width = 10, height = 8.9)
    plot(log10(df$baseMean), df$log2FoldChange, pch = 20,
         col = ifelse(is_sig, "firebrick", "grey70"),
         xlab = "log10 baseMean", ylab = "log2 fold change",
         main = paste("MA:", gsub("_", " vs ", nm)))
    abline(h = 0, col = "grey30")
    dev.off()
  }

  for (comp in comparisons) plot_one(comp$df, comp$label)

  # run metadata file for reproducibility
  write.csv(data.frame(
    formula = formula,
    contrasts = paste(vapply(contrasts, function(x) {
      paste(x$numerator, "vs", x$denominator)
    }, ""), collapse = "; "),
    fdr = fdr,
    log2fc = log2fc,
    min_count = min_count,
    min_samples = min_samples,
    genes_file = if (is.null(gene_list_path)) "" else gene_list_path
  ), file.path(outdir, "deseq2_run_meta.csv"), row.names = FALSE)

  cat("\n=== DESeq2 analysis complete! ===\n")
  cat("Output directory:", outdir, "\n")
  invisible(list(dds = dds, normalized_counts = norm_counts, results = comparisons))
}

# ------------------------------------------------------------------
# CLI wrapper
# ------------------------------------------------------------------
if (sys.nframe() == 0 && !interactive()) {
  opts <- parse_args()
  run_standard_deseq2(
    count_matrix_path = opts$count_matrix,
    meta_path = opts$meta,
    formula = opts$formula,
    contrasts = opts$contrasts,
    outdir = opts$outdir,
    fdr = opts$fdr,
    log2fc = opts$log2fc,
    min_count = opts$min_count,
    min_samples = opts$min_samples,
    gene_list_path = opts$genes,
    alpha = opts$alpha,
    shrink = opts$shrink
  )
}
