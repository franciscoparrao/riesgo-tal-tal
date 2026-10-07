# Robustez III (review ciego Natural Hazards):
#   1. Índice de riesgo sin escala: percentil del índice de alcance entre las unidades × consecuencias
#   2. Bootstrap espacial (bloques de 1 km) del AUC de HAND en cada comuna
#   3. Incertidumbre de la selección del modelo de alcance (LOCO): bootstrap de manzanas dentro de
#      cada comuna, diferencia HAND vs HAND + distancia y frecuencia con que gana cada uno
#   4. Taltal separado: manzanas que drenan a un cauce vs manzanas costeras (HAND = elevación)
#   5. Corrección de Holm para las diferencias pareadas (Tabla S3)

riesgo_sin_escala <- function(d, semilla = 42) {
  pct <- (rank(d$p_alcance) - 0.5) / nrow(d)
  C <- (d$v_extendida_indice + d$exposicion + d$locc) / 3
  capas <- list("Reach index (logistic)" = d$p_alcance, "Reach percentile" = pct,
                "Risk, logistic reach x C" = d$p_alcance * C, "Risk, reach percentile x C" = pct * C)
  sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))
  do.call(rbind, lapply(names(capas), function(n) {
    r <- auc_ic(capas[[n]], d$afectada, semilla = semilla)
    data.frame(capa = n, auc = round(r[["auc"]], 3), ic95_inf = round(r[["ic_inf"]], 3),
               ic95_sup = round(r[["ic_sup"]], 3), rho_pob_afectada = round(sp(capas[[n]], d$pob_afectada), 3),
               rho_con_alcance = round(sp(capas[[n]], d$p_alcance), 3))
  }))
}

bloques_con_coordenadas <- function(catastro, censo_gpkg, area_2015, p, umbral = 2000) {
  do.call(rbind, lapply(seq_len(nrow(comunas_validacion)), function(i) {
    comuna <- comunas_validacion$comuna[i]
    m <- manzanas_con_distancia(comuna, comunas_validacion$cut[i], catastro, censo_gpkg, area_2015, p)
    if (comuna != "TALTAL") m <- m[m$dist_area <= 5000, ]
    dir <- dirname(salida(p, "validacion_externa", carpeta_ascii(comuna), "x"))
    xy <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(m)))
    data.frame(comuna = comuna, x = xy[, 1], y = xy[, 2], afectada = m$afectada,
               hand = extraer(m, file.path(dir, sprintf("hand_costa_%d.tif", umbral))),
               hand_cruda = extraer(m, file.path(dir, sprintf("hand_%d.tif", umbral))),
               dist = extraer(m, file.path(dir, sprintf("dist_cauce_%d.tif", umbral))))
  }))
}

auc_bloques_comunas <- function(b, semilla = 42, bloque_m = 1000) {
  do.call(rbind, lapply(split(b, b$comuna), function(x) {
    bl <- paste(floor(x$x / bloque_m), floor(x$y / bloque_m))
    h <- auc_bloques(-x$hand, x$afectada, bl, semilla = semilla)
    dd <- auc_bloques(-x$dist, x$afectada, bl, semilla = semilla)
    data.frame(comuna = x$comuna[1], n = nrow(x), afectadas = sum(x$afectada), bloques_1km = h[["n_bloques"]],
               auc_hand = round(h[["auc"]], 3), ic_hand_bloques = paste(round(h[["ic_inf"]], 2), round(h[["ic_sup"]], 2), sep = "-"),
               auc_dist = round(dd[["auc"]], 3), ic_dist_bloques = paste(round(dd[["ic_inf"]], 2), round(dd[["ic_sup"]], 2), sep = "-"))
  }))
}

seleccion_bootstrap <- function(d, calibracion, B = 200, semilla = 42) {
  dc <- d[d$comuna %in% calibracion & is.finite(d$hand) & is.finite(d$dist), ]
  loco_media <- function(x, f) {
    x$w <- as.numeric(1 / table(x$comuna)[x$comuna]); x$w <- x$w / mean(x$w)
    mean(vapply(calibracion, function(cm) {
      tr <- x[x$comuna != cm, ]; te <- x[x$comuna == cm, ]
      fit <- suppressWarnings(stats::glm(f, stats::binomial(), tr, weights = w))
      auc(stats::predict(fit, te), te$afectada)
    }, 0))
  }
  f_h <- afectada ~ log1p(hand); f_hd <- afectada ~ log1p(hand) + log1p(dist); f_d <- afectada ~ log1p(dist)
  set.seed(semilla)
  por <- split(seq_len(nrow(dc)), dc$comuna)
  res <- t(replicate(B, {
    i <- unlist(lapply(por, function(ix) sample(ix, replace = TRUE)), use.names = FALSE)
    x <- dc[i, ]
    c(hand = loco_media(x, f_h), hand_dist = loco_media(x, f_hd), dist = loco_media(x, f_d))
  }))
  res <- res[stats::complete.cases(res), , drop = FALSE]   # réplicas sin afectadas en alguna comuna (Alto del Carmen)
  obs <- c(hand = loco_media(dc, f_h), hand_dist = loco_media(dc, f_hd), dist = loco_media(dc, f_d))
  gana <- table(factor(colnames(res)[apply(res, 1, which.max)], colnames(res)))
  q <- function(v) round(stats::quantile(v, c(0.025, 0.975)), 3)
  data.frame(modelo = names(obs), loco_media = round(obs, 3),
             ic95 = vapply(colnames(res), function(k) paste(q(res[, k]), collapse = "-"), ""),
             dif_vs_hand_ic95 = vapply(colnames(res), function(k) paste(q(res[, "hand"] - res[, k]), collapse = "-"), ""),
             pct_seleccionado = round(100 * as.numeric(gana[names(obs)]) / nrow(res), 1), replicas_validas = nrow(res))
}

taltal_por_tipo <- function(b, semilla = 42) {
  x <- b[b$comuna == "TALTAL", ]
  costa <- !is.finite(x$hand_cruda)
  f <- function(v, cl) { r <- auc_ic(v, cl, semilla = semilla); sprintf("%.2f (%.2f-%.2f)", r[["auc"]], r[["ic_inf"]], r[["ic_sup"]]) }
  data.frame(grupo = c("drain to a channel (HAND)", "drain to the sea (HAND = elevation)", "all"),
             n = c(sum(!costa), sum(costa), nrow(x)), afectadas = c(sum(x$afectada[!costa]), sum(x$afectada[costa]), sum(x$afectada)),
             auc_hand = c(f(-x$hand[!costa], x$afectada[!costa]), f(-x$hand[costa], x$afectada[costa]), f(-x$hand, x$afectada)),
             auc_dist = c(f(-x$dist[!costa], x$afectada[!costa]), f(-x$dist[costa], x$afectada[costa]), f(-x$dist, x$afectada)))
}

holm_pareado <- function(csv) {
  t <- utils::read.csv(csv)
  p2 <- pmin(1, 2 * pmin(t$p_unilateral_bloque, 1 - t$p_unilateral_bloque))
  t$p_bilateral_bloque <- round(p2, 3)
  t$p_holm <- round(stats::p.adjust(p2, "holm"), 3)
  t
}
