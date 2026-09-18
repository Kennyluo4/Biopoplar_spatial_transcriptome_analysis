# Portable paths. Never point this at the original analysis directory for a rerun:
# scripts write results below POPLAR_ATLAS_ROOT using the documented study layout.
ATLAS_CODE_ROOT <- normalizePath(.code_root, mustWork = TRUE)
ATLAS_WORK_ROOT <- Sys.getenv("POPLAR_ATLAS_ROOT", file.path(ATLAS_CODE_ROOT, "work"))

atlas_input <- function(relative) {
  p <- file.path(ATLAS_WORK_ROOT, relative)
  if (!file.exists(p)) stop("Missing input: ", relative,
    ". Populate the study work directory; see docs/INPUTS.md.", call. = FALSE)
  p
}

atlas_annotations <- function(object, tissue) {
  f <- file.path(ATLAS_CODE_ROOT, "resources", "annotations", paste0(tissue, "_spot_annotations.csv"))
  a <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (anyDuplicated(a$spot)) stop("Duplicate spot IDs in curated annotation table")
  if (!all(a$spot %in% colnames(object))) stop("Rebuilt spot IDs do not match the deposited curated annotation. Check section splitting and QC before continuing.")
  object <- subset(object, cells = a$spot)
  idx <- match(colnames(object), a$spot)
  for (nm in setdiff(names(a), "spot")) object[[nm]] <- a[[nm]][idx]
  Seurat::Idents(object) <- object$celltypes
  object
}

atlas_dispatch <- function(stages) {
  args <- commandArgs(trailingOnly = TRUE)
  stage <- if (length(args)) args[1] else "--list"
  if (stage %in% c("--list", "--help")) {
    cat("Available stages:\n", paste(names(stages), collapse = "\n"),
        "\nRun: Rscript <tissue_script.R> <stage>\nSee README.md for input requirements and execution order.\n", sep = "")
    return(invisible(NULL))
  }
  if (!stage %in% names(stages)) stop("Unknown stage: ", stage)
  if (!dir.exists(ATLAS_WORK_ROOT)) stop("Create and populate POPLAR_ATLAS_ROOT first; see docs/INPUTS.md")
  ATLAS_WORK_ROOT <<- normalizePath(ATLAS_WORK_ROOT, mustWork = TRUE)
  old <- setwd(ATLAS_WORK_ROOT)
  on.exit(setwd(old), add = TRUE)
  source(file.path(ATLAS_CODE_ROOT, "R", "common_helpers.R"), local = .GlobalEnv)
  for (d in c("QC", "section_split", "clustering", "marker", "saved_obj", "tables", "trajectory", "go", "nmf", "figures", "logs")) dir.create(d, recursive=TRUE, showWarnings=FALSE)
  set.seed(12345)
  on.exit(writeLines(capture.output(sessionInfo()), file.path("logs", paste0(basename(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])), "_",stage,"_sessionInfo.txt"))), add=TRUE)
  stages[[stage]]()
  invisible(NULL)
}
