# Load and install all libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", ask = FALSE, quiet = TRUE)
BiocManager::install("edgeR", ask = FALSE, quiet = TRUE)
BiocManager::install("apeglm", ask = FALSE, quiet = TRUE)
BiocManager::install('EnhancedVolcano', ask = FALSE, quiet = TRUE)
BiocManager::install("pasilla", ask = FALSE, quiet = TRUE)
BiocManager::install("DEGreport", ask = FALSE, quiet = TRUE)

if (!require("factoextra", quietly = TRUE))
  install.packages("factoextra")

if (!require("purrr", quietly = TRUE))
  install.packages("purrr")

if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

library(DESeq2)
library(edgeR)
library(readr)
library(factoextra)
library(purrr)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(styler)
library(knitr)
library(DEGreport)
library(ggplot2)
library(genefilter)
library(grid)