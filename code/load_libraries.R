# Load all libraries required by the report.
#
# Package installation is intentionally kept out of render/preview. Install the
# conda environment from environment.yml, then render the document from that
# environment. This keeps Quarto from silently compiling packages during setup.
message("Loading libraries")

attached_packages <- c(
  "DESeq2",
  "edgeR",
  "readr",
  "dplyr",
  "knitr",
  "ggplot2",
  "grid"
)

namespace_packages <- c(
  attached_packages,
  "apeglm",
  "circlize",
  "ComplexHeatmap",
  "GenomicRanges",
  "IRanges",
  "limma",
  "S4Vectors",
  "scales"
)

optional_namespace_packages <- c("rtracklayer")

missing_packages <- namespace_packages[
  !vapply(namespace_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ".\nCreate/update the conda environment first:\n",
      "  conda env update -f environment.yml --prune\n",
      "  conda activate differential-interaction-analysis"
    ),
    call. = FALSE
  )
}

invisible(lapply(attached_packages, library, character.only = TRUE))

missing_optional_packages <- optional_namespace_packages[
  !vapply(optional_namespace_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_optional_packages) > 0) {
  message(
    "Optional R packages not available: ",
    paste(missing_optional_packages, collapse = ", "),
    ". GFF annotation will use the built-in parser fallback where possible."
  )
}
