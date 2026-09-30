#!/usr/bin/env Rscript
# Standardized WGCNA wrapper for UGT co-expression mining
#
# Required inputs:
#   --expr     gene-level expression matrix (gene_id in first column, one sample per column)
#   --meta     sample metadata TSV: at least two columns named sample and trait
#   --bait     bait gene id present in the expression matrix
#   --cands    candidate gene list, one gene id per line (optional)
#
# Optional:
#   --outdir     output directory (default: wgcna_results)
#   --network    network type, unsigned/signed/signed hybrid (default: unsigned)
#   --power      fixed soft power (default: auto)
#   --min-size   minimum module size (default: 30)
#   --merge       module merge cut height (default: 0.25)
#   --max-block   max block size for blockwiseModules (default: 30000)
#   --threshold   expression filter threshold (default: 1)
#   --min-samples minimum number of samples above threshold (default: NULL)
#   --nthreads    WGCNA threads (default: 10)
#   --save-tom    save TOM file, 0/1 (default: 0)
#
# Example:
#   Rscript run_wgcna.R \
#     --expr  gene_TPM_matrix.tsv \
#     --meta  sample_meta.tsv \
#     --bait  MSTRG.7391 \
#     --cands ugt_candidates.txt \
#     --outdir results/wgcna

suppressPackageStartupMessages({
  library(WGCNA)
  library(pheatmap)
  library(fastcluster)
})

# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- list(
    expr = NULL, meta = NULL, bait = NULL, cands = NULL,
    outdir = "wgcna_results",
    network = "unsigned",
    power = NA,
    min_size = 30,
    merge = 0.25,
    max_block = 30000,
    threshold = 1.0,
    min_samples = NULL,
    nthreads = 10,
    save_tom = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    take_val <- key %in% c(
      "--expr", "--meta", "--bait", "--cands", "--outdir", "--network",
      "--power", "--min-size", "--merge", "--max-block", "--threshold",
      "--min-samples", "--nthreads", "--save-tom"
    )
    if (!take_val) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    switch(key,
      "--expr"        = { opts$expr <- val; i <- i + 2 },
      "--meta"        = { opts$meta <- val; i <- i + 2 },
      "--bait"        = { opts$bait <- val; i <- i + 2 },
      "--cands"       = { opts$cands <- val; i <- i + 2 },
      "--outdir"      = { opts$outdir <- val; i <- i + 2 },
      "--network"     = { opts$network <- val; i <- i + 2 },
      "--power"       = { opts$power <- suppressWarnings(as.numeric(val)); i <- i + 2 },
      "--min-size"    = { opts$min_size <- as.integer(val); i <- i + 2 },
      "--merge"       = { opts$merge <- as.numeric(val); i <- i + 2 },
      "--max-block"   = { opts$max_block <- as.integer(val); i <- i + 2 },
      "--threshold"   = { opts$threshold <- as.numeric(val); i <- i + 2 },
      "--min-samples" = { opts$min_samples <- as.integer(val); i <- i + 2 },
      "--nthreads"    = { opts$nthreads <- as.integer(val); i <- i + 2 },
      "--save-tom"    = { opts$save_tom <- as.logical(as.integer(val)); i <- i + 2 }
    )
  }

  required <- c("expr", "meta", "bait")
  missing <- required[sapply(required, function(x) is.null(opts[[x]]) || is.na(opts[[x]]))]
  if (length(missing) > 0) {
    stop(
      "Missing required argument(s): ", paste0("--", missing, collapse = ", "),
      "\nRun: Rscript run_wgcna.R --help",
      call. = FALSE
    )
  }
  opts
}

# ------------------------------------------------------------------
# Input loading helpers
# ------------------------------------------------------------------
read_expr <- function(path) {
  if (!file.exists(path)) {
    stop("Expression matrix not found: ", path, call. = FALSE)
  }
  x <- read.delim(path, header = TRUE, row.names = 1, check.names = FALSE, quote = "")
  num_cols <- sapply(x, is.numeric)
  if (!all(num_cols)) {
    stop(
      "Expression matrix columns must all be numeric; non-numeric columns: ",
      paste(names(num_cols)[!num_cols], collapse = ", "),
      call. = FALSE
    )
  }
  x
}

read_metadata <- function(path) {
  if (!file.exists(path)) {
    stop("Sample metadata not found: ", path, call. = FALSE)
  }
  m <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% colnames(m)) {
    stop("Metadata must contain a 'sample' column.", call. = FALSE)
  }
  if (!"trait" %in% colnames(m)) {
    stop("Metadata must contain a 'trait' column.", call. = FALSE)
  }
  rownames(m) <- m$sample
  m$trait <- factor(m$trait)
  m
}

read_candidate_list <- function(path) {
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
run_ugt_wgcna <- function(expr_path, meta_path, bait_gene,
                          cand_path = NULL,
                          outdir = "wgcna_results",
                          network_type = "unsigned",
                          soft_power = NA,
                          min_module_size = 30,
                          merge_cut_height = 0.25,
                          max_block_size = 30000,
                          expr_threshold = 1,
                          min_expr_samples = NULL,
                          nthreads = 10,
                          save_tom = FALSE) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  enableWGCNAThreads(nthreads)

  exprData <- read_expr(expr_path)
  meta <- read_metadata(meta_path)
  ugtIDs <- read_candidate_list(cand_path)

  common_samples <- intersect(colnames(exprData), rownames(meta))
  if (length(common_samples) == 0) {
    stop("No samples in expression matrix match metadata.", call. = FALSE)
  }
  exprData <- exprData[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  meta <- meta[colnames(exprData), , drop = FALSE]

  cat("Expression matrix:", nrow(exprData), "genes x", ncol(exprData), "samples\n")
  cat("Samples used:", paste(colnames(exprData), collapse = ", "), "\n")
  cat("Candidate genes listed:", length(ugtIDs), "\n")

  if (is.null(min_expr_samples)) {
    min_expr_samples <- max(5, round(ncol(exprData) / 2))
  }
  keep <- apply(exprData, 1, function(x) sum(x > expr_threshold) >= min_expr_samples)
  exprFiltered <- exprData[keep, , drop = FALSE]
  cat("Genes after expression filter (TPM >", expr_threshold,
      "in >=", min_expr_samples, "samples):", nrow(exprFiltered), "\n")

  if (!(bait_gene %in% rownames(exprFiltered))) {
    if (bait_gene %in% rownames(exprData)) {
      exprFiltered <- rbind(exprFiltered, exprData[bait_gene, , drop = FALSE])
      cat("Bait gene filtered out normally; added back for network.\n")
    } else {
      stop("Bait gene '", bait_gene, "' not found in expression matrix.", call. = FALSE)
    }
  }

  datExpr0 <- t(log2(exprFiltered + 1))
  cat("WGCNA expression matrix:", nrow(datExpr0), "samples x",
      ncol(datExpr0), "genes\n")

  # binary trait encoding for numerical module-trait correlation
  trait_levels <- levels(meta$trait)
  traits_binary <- sapply(trait_levels, function(lv) as.numeric(meta$trait == lv))
  if (!is.matrix(traits_binary)) {
    traits_binary <- matrix(traits_binary, ncol = 1)
  }
  colnames(traits_binary) <- trait_levels
  rownames(traits_binary) <- rownames(meta)

  cat("Traits:\n")
  print(traits_binary)

  # sample QC
  sampleTree <- hclust(dist(datExpr0), method = "average")
  pdf(file.path(outdir, "sample_clustering.pdf"), width = 12, height = 8)
  plot(sampleTree, main = "Sample clustering", xlab = "", sub = "",
       cex = 0.8, hang = -1)
  traitColors <- numbers2colors(traits_binary, signed = FALSE)
  plotDendroAndColors(sampleTree, traitColors,
                      groupLabels = colnames(traits_binary),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()

  # soft threshold
  powers_used <- 1:20
  sft <- pickSoftThreshold(datExpr0, powerVector = powers_used,
                           networkType = network_type, verbose = 0)

  pdf(file.path(outdir, "soft_threshold.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale Independence")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers_used, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = "Mean Connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = powers_used, cex = cex1, col = "red")
  dev.off()

  if (is.na(soft_power)) {
    soft_power <- sft$powerEstimate
    if (is.na(soft_power)) {
      cat("Auto soft threshold is NA, using 6\n")
      soft_power <- 6
    }
  }
  cat("Selected soft power:", soft_power, "\n")

  net <- blockwiseModules(
    datExpr0,
    power = soft_power,
    networkType = network_type,
    TOMType = network_type,
    minModuleSize = min(min_module_size, ncol(datExpr0)),
    mergeCutHeight = merge_cut_height,
    numericLabels = TRUE,
    maxBlockSize = max_block_size,
    saveTOMs = save_tom,
    saveTOMFileBase = file.path(outdir, "TOM"),
    verbose = 0
  )

  moduleColors <- labels2colors(net$colors)
  geneNames <- colnames(datExpr0)

  cat("Modules identified:", length(unique(net$colors)), "\n")
  print(table(moduleColors))
  if (length(unique(moduleColors)) == 1 && unique(moduleColors) == "grey") {
    cat("NOTE: no non-grey modules were detected. ",
        "This usually means the dataset is too small or too many genes were filtered. ",
        "Try a larger expression matrix or lower --min-size.\n")
  }

  # module dendrogram (all blocks, if single block use net$dendrograms[[1]])
  pdf(file.path(outdir, "module_dendrogram.pdf"), width = 12, height = 8)
  if (!is.null(net$dendrograms) && length(net$dendrograms) >= 1 &&
      !is.null(net$blockGenes) && length(net$blockGenes) >= 1) {
    first_block_genes <- net$blockGenes[[1]]
    first_block_labels <- first_block_genes[first_block_genes > 0]
    if (length(first_block_labels) > 0) {
      plotDendroAndColors(
        net$dendrograms[[1]],
        moduleColors[first_block_labels],
        "Module colors",
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05
      )
    }
  }
  dev.off()

  # module eigengenes and module-trait correlation
  MEs0 <- moduleEigengenes(datExpr0, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, traits_binary, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr0))

  mt_df <- data.frame(
    Module = rownames(moduleTraitCor),
    round(moduleTraitCor, 4),
    round(moduleTraitPvalue, 4),
    stringsAsFactors = FALSE
  )
  colnames(mt_df)[(ncol(traits_binary) + 2):ncol(mt_df)] <-
    paste0(colnames(traits_binary), "_pval")
  write.csv(mt_df, file.path(outdir, "module_trait_correlation.csv"),
            row.names = FALSE)

  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)
  pdf(file.path(outdir, "module_trait_heatmap.pdf"), width = 10, height = 12)
  par(mar = c(6, 10, 3, 3))
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(traits_binary),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    zlim = c(-1, 1),
    main = "Module-trait relationships"
  )
  dev.off()

  # locate bait gene
  baitIdx <- which(geneNames == bait_gene)
  if (length(baitIdx) == 0) {
    partial <- grep(bait_gene, geneNames, value = TRUE)
    if (length(partial) > 0) {
      bait_gene <- partial[1]
      baitIdx <- which(geneNames == bait_gene)
      cat("Found partial bait match:", bait_gene, "\n")
    }
  }

  if (length(baitIdx) == 0) {
    cat("FATAL: Could not find bait gene in expression data.\n")
  } else {
    baitModule <- moduleColors[baitIdx]
    cat("Bait gene:", bait_gene, "\n")
    cat("Bait gene module:", baitModule, "\n")

    coExprGenes <- geneNames[moduleColors == baitModule]
    cat("Co-expression partners:", length(coExprGenes), "genes\n")
    writeLines(sort(coExprGenes), file.path(outdir, "coexpression_partners.txt"))

    baitME <- paste0("ME", baitModule)
    if (baitME %in% colnames(MEs)) {
      cat("Bait module trait correlations:\n")
      for (t in colnames(traits_binary)) {
        r <- cor(MEs[[baitME]], traits_binary[, t], use = "p")
        p <- corPvalueStudent(r, nSamples = nrow(datExpr0))
        cat(sprintf("  %s: r=%.3f, p=%.4g\n", t, r, p))
      }
    }

    ugtInModule <- intersect(ugtIDs, coExprGenes)
    cat("UGT candidates in bait module:", length(ugtInModule), "\n")
    if (length(ugtInModule) > 0) {
      writeLines(sort(ugtInModule), file.path(outdir, "ugt_in_bait_module.txt"))
      cat("UGTs in module:\n")
      print(sort(ugtInModule))
    }

    # Pearson of bait-module UGTs vs bait
    exprLog <- log2(exprFiltered + 1)
    if (bait_gene %in% rownames(exprLog)) {
      baitVec <- as.numeric(exprLog[bait_gene, ])
      cand <- setdiff(ugtInModule, bait_gene)
      cand <- cand[cand %in% rownames(exprLog)]
      pr <- data.frame(
        gene = character(0), module = character(0),
        r_to_bait = numeric(0), pval = character(0),
        stringsAsFactors = FALSE
      )
      for (g in cand) {
        r <- cor(baitVec, as.numeric(exprLog[g, ]), method = "pearson")
        p <- corPvalueStudent(r, nSamples = length(baitVec))
        pr <- rbind(pr, data.frame(
          gene = g, module = baitModule, r_to_bait = round(r, 4),
          pval = format.pval(p, digits = 3), stringsAsFactors = FALSE
        ))
      }
      if (nrow(pr) > 0) pr <- pr[order(-pr$r_to_bait), , drop = FALSE]
      write.csv(pr, file.path(outdir, "candidate_ugt_pearson.csv"),
                row.names = FALSE)
      cat("Candidates ranked by r_to_bait:", nrow(pr),
          "(bait module UGTs, bait excluded)\n")
    }

    if (length(ugtIDs) > 0) {
      ugtInExpr <- ugtIDs[ugtIDs %in% geneNames]
      ugtModuleMap <- data.frame(
        UGT = ugtInExpr,
        Module = moduleColors[match(ugtInExpr, geneNames)],
        stringsAsFactors = FALSE
      )
      ugtModuleMap <- ugtModuleMap[order(ugtModuleMap$Module), , drop = FALSE]
      write.csv(ugtModuleMap, file.path(outdir, "ugt_to_module_map.csv"),
                row.names = FALSE)
      cat("UGT-module distribution:\n")
      print(table(ugtModuleMap$Module))

      ugtInFiltered <- ugtIDs[ugtIDs %in% rownames(exprFiltered)]
      if (length(ugtInFiltered) > 0) {
        ugtTPM <- exprFiltered[ugtInFiltered, , drop = FALSE]
        write.csv(ugtTPM, file.path(outdir, "ugt_tpm_matrix.csv"))
      }
    }

    if (length(ugtInModule) > 0) {
      ugtModuleExpr <- exprFiltered[ugtInModule, , drop = FALSE]
      ugtModuleExpr <- log2(ugtModuleExpr + 1)
      if (!(bait_gene %in% rownames(ugtModuleExpr)) &&
          bait_gene %in% rownames(exprFiltered)) {
        ugtModuleExpr <- rbind(ugtModuleExpr,
                               log2(exprFiltered[bait_gene, , drop = FALSE] + 1))
      }
      pdf(file.path(outdir, "ugt_bait_module_heatmap.pdf"), width = 10, height = 8)
      pheatmap(
        ugtModuleExpr,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "average",
        show_rownames = TRUE,
        fontsize_row = 8,
        main = paste0("UGT candidates in ", baitModule, " module")
      )
      dev.off()
    }

    if (baitME %in% colnames(MEs)) {
      kME <- cor(datExpr0, MEs[[baitME]], use = "p")
      moduleGenes <- which(moduleColors == baitModule)
      hubGenes <- data.frame(
        Gene = geneNames[moduleGenes],
        kME = kME[moduleGenes, 1],
        isUGT = geneNames[moduleGenes] %in% ugtIDs,
        stringsAsFactors = FALSE
      )
      hubGenes <- hubGenes[order(-hubGenes$kME), , drop = FALSE]
      write.csv(hubGenes, file.path(outdir, "bait_module_hub_genes.csv"),
                row.names = FALSE)
      cat("\nTop 20 hub genes in bait module:\n")
      print(head(hubGenes, 20))
    }
  }

  moduleAssign <- data.frame(
    Gene = geneNames,
    Module = moduleColors,
    stringsAsFactors = FALSE
  )
  write.csv(moduleAssign, file.path(outdir, "gene_module_assignments.csv"),
            row.names = FALSE)

  # simple summary row for convenience
  write.csv(data.frame(
    soft_power = soft_power,
    genes_used = nrow(exprFiltered),
    samples_used = ncol(exprFiltered),
    modules = length(unique(net$colors)),
    bait_gene = bait_gene,
    bait_module = if (length(baitIdx) > 0) moduleColors[baitIdx] else NA,
    bait_module_ugts = if (exists("ugtInModule") && !is.null(ugtInModule)) length(ugtInModule) else NA
  ), file.path(outdir, "wgcna_summary.csv"), row.names = FALSE)

  cat("\n=== WGCNA analysis complete! ===\n")
  cat("Output directory:", outdir, "\n")
  invisible(list(
    expr_filtered = exprFiltered,
    datExpr = datExpr0,
    module_colors = moduleColors,
    module_eigengenes = MEs,
    module_trait_cor = moduleTraitCor,
    module_trait_pval = moduleTraitPvalue,
    soft_power = soft_power
  ))
}

# ------------------------------------------------------------------
# CLI wrapper
# ------------------------------------------------------------------
if (sys.nframe() == 0 && !interactive()) {
  opts <- parse_args()
  run_ugt_wgcna(
    expr_path = opts$expr,
    meta_path = opts$meta,
    bait_gene = opts$bait,
    cand_path = opts$cands,
    outdir = opts$outdir,
    network_type = opts$network,
    soft_power = opts$power,
    min_module_size = opts$min_size,
    merge_cut_height = opts$merge,
    max_block_size = opts$max_block,
    expr_threshold = opts$threshold,
    min_expr_samples = opts$min_samples,
    nthreads = opts$nthreads,
    save_tom = opts$save_tom
  )
}
