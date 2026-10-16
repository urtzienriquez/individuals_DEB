#!/usr/bin/env Rscript
# merge chunk RDS files into final output

args <- commandArgs(trailingOnly = TRUE)
SCENARIO_FILTER <- NULL

for (a in args) {
  if (grepl("^--scenario=", a)) {
    SCENARIO_FILTER <- sub("^--scenario=", "", a)
  }
}

results_dir <- "../results"

scenarios <- if (!is.null(SCENARIO_FILTER)) {
  SCENARIO_FILTER
} else {
  c("baseline", "w2", "w4")
}

for (sc in scenarios) {
  cat(sprintf("\nMerging scenario: %s\n", sc))

  pattern <- sprintf("europe_results_%s_chunk.*of.*\\.rds$", sc)
  files <- list.files(results_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    cat(sprintf("  No chunk files found matching '%s' in %s - skipping\n", pattern, results_dir))
    next
  }

  cat(sprintf("  Found %d chunk files\n", length(files)))

  chunk_nums <- sort(unique(as.integer(
    sub(sprintf(".*chunk(\\d+)of\\d+\\.rds$", sc), "\\1", basename(files))
  )))
  n_chunks <- max(as.integer(
    sub(sprintf(".*chunk\\d+of(\\d+)\\.rds$", sc), "\\1", basename(files))
  ))
  cat(sprintf("  Chunks present: %s (of %d expected)\n",
              paste(chunk_nums, collapse = ", "), n_chunks))

  missing <- setdiff(seq_len(n_chunks), chunk_nums)
  if (length(missing) > 0) {
    cat(sprintf("  WARNING: missing chunks %s\n", paste(missing, collapse = ", ")))
  }

  cat("  Reading chunks ... ")
  chunks <- lapply(files, readRDS)
  results <- do.call(rbind, chunks)
  rownames(results) <- NULL
  cat("done\n")

  n_ok <- sum(results$type != "error" | is.na(results$type))
  cat(sprintf("  Total rows: %d\n", nrow(results)))
  cat(sprintf("  Successful cells: %d\n", n_ok))

  out_file <- file.path(results_dir, sprintf("europe_results_%s.rds", sc))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(results, out_file)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\nDone.\n")
