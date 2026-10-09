#!/usr/bin/env Rscript
# Controles previos al envío (review ciego NH):
#  1. HAND de SurtGIS vs. WhiteboxTools 2.4 (bin/hand_whitebox.py) en las manzanas de las seis
#     comunas: correlación de Spearman, AUC por comuna, LOCO y AUC del modelo transferido a Taltal.
#  2. GFI con exponente n fijo (no ajustado): LOCO en Atacama y AUC en Taltal para n = 0-0.6.
# Uso (desde la raíz, con el pipeline corrido):  Rscript bin/control_hand_gfi.R
suppressPackageStartupMessages({ library(terra) })
Sys.setenv(RIESGO_MODO = "corregido")
targets::tar_source("R")
st <- "_targets/corregido"
b <- targets::tar_read(bloques_amp, store = st)
b <- b[stats::complete.cases(b[, c("hand", "area_cauce")]), ]
dir_ve <- "salidas/corregido/validacion_externa"
out <- "salidas/corregido/robustez"

# --- 1. HAND Whitebox con la misma regla costera (NA -> elevación) ---
b$hand_wbt <- NA_real_
for (cm in unique(b$comuna)) {
  i <- b$comuna == cm
  d <- file.path(dir_ve, carpeta_ascii(cm), "wbt")
  h <- rast(file.path(d, "hand_wbt_2000.tif")); z <- rast(file.path(d, "dem_fill.tif"))
  h <- cover(h, clamp(z, lower = 0))
  b$hand_wbt[i] <- extract(h, cbind(b$x[i], b$y[i]))[, 1]
}
b <- b[!is.na(b$hand_wbt), ]
cal <- setdiff(unique(b$comuna), "TALTAL")
por_comuna <- do.call(rbind, lapply(split(b, b$comuna), function(x) data.frame(
  comuna = x$comuna[1], n = nrow(x), rho_spearman = round(stats::cor(x$hand, x$hand_wbt, method = "spearman"), 3),
  dif_mediana_m = round(stats::median(x$hand_wbt - x$hand), 2),
  auc_surtgis = round(auc(-x$hand, x$afectada), 3), auc_whitebox = round(auc(-x$hand_wbt, x$afectada), 3))))
# LOCO y transferencia con el mismo modelo logístico (comunas con igual peso)
ajuste <- function(tr, v) {
  tr$w <- as.numeric(1 / table(tr$comuna)[tr$comuna]); tr$w <- tr$w / mean(tr$w)
  tr$h <- tr[[v]]
  suppressWarnings(stats::glm(afectada ~ log1p(h), stats::binomial(), tr, weights = w))
}
loco <- function(v) mean(vapply(cal, function(cm) {
  te <- b[b$comuna == cm, ]; te$h <- te[[v]]
  auc(stats::predict(ajuste(b[b$comuna %in% setdiff(cal, cm), ], v), te), te$afectada)
}, 0))
tt <- b[b$comuna == "TALTAL", ]
transf <- function(v) { te <- tt; te$h <- te[[v]]; auc(stats::predict(ajuste(b[b$comuna %in% cal, ], v), te), te$afectada) }
resumen_hand <- data.frame(
  software = c("SurtGIS 1.5.1", "WhiteboxTools 2.4.0"),
  rho_global = c(1, round(stats::cor(b$hand, b$hand_wbt, method = "spearman"), 3)),
  loco_media = round(c(loco("hand"), loco("hand_wbt")), 3),
  auc_taltal_bloques = round(c(transf("hand"), transf("hand_wbt")), 3),
  beta1 = round(c(stats::coef(ajuste(b[b$comuna %in% cal, ], "hand"))[2], stats::coef(ajuste(b[b$comuna %in% cal, ], "hand_wbt"))[2]), 3))
write.csv(por_comuna, file.path(out, "hand_whitebox_por_comuna.csv"), row.names = FALSE)
write.csv(resumen_hand, file.path(out, "hand_whitebox_resumen.csv"), row.names = FALSE)

# --- 2. GFI con n fijo (sin ajuste; el puntaje se usa directo) ---
gfi_fijo <- do.call(rbind, lapply(seq(0, 0.6, 0.05), function(n) {
  s <- gfi(b, n)
  a <- vapply(cal, function(cm) { i <- b$comuna == cm; auc(s[i], b$afectada[i]) }, 0)
  i <- b$comuna == "TALTAL"
  data.frame(n = n, loco_media = round(mean(a), 3), loco_min = round(min(a), 3),
             auc_taltal = round(auc(s[i], b$afectada[i]), 3))
}))
write.csv(gfi_fijo, file.path(out, "gfi_n_fijo.csv"), row.names = FALSE)

print(por_comuna); print(resumen_hand); print(gfi_fijo)
