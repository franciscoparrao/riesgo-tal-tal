#!/usr/bin/env Rscript
# Uso: Rscript run.R [corregido]
# Store de targets en _targets/<modo>/ y resultados en salidas/<modo>/.

modo <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(modo)) modo <- "corregido"
Sys.setenv(RIESGO_MODO = modo)
targets::tar_make(store = file.path("_targets", modo))
