# Selección ampliada del término de alcance y contrastes conjuntos (review ciego NH, ronda 2):
#   1. Candidatos: HAND, distancia al cauce, HAND + distancia, elevación, pendiente y el Geomorphic
#      Flood Index (GFI = n ln A_r - ln(1 + HAND); Samela et al. 2017), con n elegido en las comunas
#      de entrenamiento. Selección por AUC dejando una comuna fuera (LOCO) en Atacama.
#   2. Contrastes con las seis comunas juntas: AUC media dentro de comunas (pesos iguales), diferencias
#      pareadas con bootstrap de bloques de 1 km dentro de cada comuna, Holm; y modelos logísticos
#      mixtos (intercepto aleatorio por comuna) comparados por AIC.
#   3. Barrido del peso de la amenaza en la agregación aditiva (Taltal).
# A_r (área aportante del cauce receptor) viene de los DEM de cuenca completa
# (salidas/corregido/movilidad/region/*/area_cauce_2000.tif; bin/area_cauce_destino.py).

region_de <- c("TALTAL" = "taltal", "CHAÑARAL" = "salado", "DIEGO DE ALMAGRO" = "salado",
               "TIERRA AMARILLA" = "copiapo", "COPIAPO" = "copiapo", "ALTO DEL CARMEN" = "huasco")

bloques_ampliado <- function(b, p, dir_region = "salidas/corregido/movilidad/region") {
  do.call(rbind, lapply(split(b, b$comuna), function(x) {
    xy <- cbind(x$x, x$y)
    dir <- dirname(salida(p, "validacion_externa", carpeta_ascii(x$comuna[1]), "x"))
    ex <- function(f) terra::extract(terra::rast(f), xy)[[1]]
    x$elev <- -ex(file.path(dir, "triviales", "elev_neg.tif"))
    x$pendiente <- -ex(file.path(dir, "triviales", "pendiente_neg.tif"))
    x$area_cauce <- ex(file.path(dir_region, region_de[[x$comuna[1]]], "area_cauce_2000.tif"))
    x
  }))
}

gfi <- function(x, n) n * log(x$area_cauce) - log1p(x$hand)

# Puntaje de cada candidato a partir de los datos de entrenamiento (tr) aplicado a te.
candidatos <- list(
  "HAND"                   = list(f = afectada ~ log1p(hand)),
  "Distance to channel"    = list(f = afectada ~ log1p(dist)),
  "HAND + distance"        = list(f = afectada ~ log1p(hand) + log1p(dist)),
  "Elevation"              = list(f = afectada ~ elev),
  "Slope"                  = list(f = afectada ~ pendiente),
  "GFI (n fitted)"         = list(gfi = TRUE)
)
auc_media_comunas <- function(score, d) mean(vapply(split(seq_len(nrow(d)), d$comuna),
                                                    function(i) auc(score[i], d$afectada[i]), 0))
puntaje <- function(cand, tr, te, ns = seq(0, 0.6, 0.05)) {
  if (isTRUE(cand$gfi)) {
    a <- vapply(ns, function(n) auc_media_comunas(gfi(tr, n), tr), 0)
    return(list(score = gfi(te, ns[which.max(a)]), n = ns[which.max(a)]))
  }
  tr$w <- as.numeric(1 / table(tr$comuna)[tr$comuna]); tr$w <- tr$w / mean(tr$w)
  fit <- suppressWarnings(stats::glm(cand$f, stats::binomial(), tr, weights = w))
  list(score = stats::predict(fit, te), n = NA)
}

seleccion_ampliada <- function(d, semilla = 42) {
  d <- d[stats::complete.cases(d[, c("hand", "dist", "elev", "pendiente", "area_cauce")]), ]
  cal <- setdiff(unique(d$comuna), "TALTAL"); dc <- d[d$comuna %in% cal, ]; tt <- d[d$comuna == "TALTAL", ]
  do.call(rbind, lapply(names(candidatos), function(nm) {
    loco <- vapply(cal, function(cm) {
      r <- puntaje(candidatos[[nm]], dc[dc$comuna != cm, ], dc[dc$comuna == cm, ])
      auc(r$score, dc$afectada[dc$comuna == cm])
    }, 0)
    rt <- puntaje(candidatos[[nm]], dc, tt)
    ic <- auc_ic(rt$score, tt$afectada, semilla = semilla)
    data.frame(candidato = nm, loco_media = round(mean(loco), 3), loco_min = round(min(loco), 3),
               loco_por_comuna = paste(sprintf("%s %.2f", substr(cal, 1, 4), loco), collapse = "; "),
               n_gfi = rt$n, auc_taltal = round(ic[["auc"]], 3),
               ic95_taltal = paste(round(ic[["ic_inf"]], 2), round(ic[["ic_sup"]], 2), sep = "-"))
  }))
}

# 2. Contrastes conjuntos: puntajes crudos (signo: mayor = más alcance), n del GFI fijado en Atacama.
contrastes_conjuntos <- function(d, B = 2000, semilla = 42, bloque_m = 1000) {
  d <- d[stats::complete.cases(d[, c("hand", "dist", "elev", "pendiente", "area_cauce")]), ]
  cal <- setdiff(unique(d$comuna), "TALTAL")
  n_gfi <- puntaje(candidatos[["GFI (n fitted)"]], d[d$comuna %in% cal, ], d)$n
  S <- data.frame(HAND = -d$hand, "Distance to channel" = -d$dist, Elevation = -d$elev,
                  Slope = -d$pendiente, GFI = gfi(d, n_gfi), check.names = FALSE)
  d$bloque <- paste(d$comuna, floor(d$x / bloque_m), floor(d$y / bloque_m))
  por <- lapply(split(seq_len(nrow(d)), d$comuna), function(i) split(i, d$bloque[i]))
  set.seed(semilla)
  rep_idx <- replicate(B, unlist(lapply(por, function(bl) unlist(bl[sample(length(bl), replace = TRUE)],
                                                                 use.names = FALSE)), use.names = FALSE),
                       simplify = FALSE)
  am <- function(s, i) { x <- d[i, ]; mean(vapply(split(seq_along(i), x$comuna), function(j)
    if (length(unique(x$afectada[j])) < 2) NA_real_ else auc(s[i][j], x$afectada[j]), 0), na.rm = TRUE) }
  obs <- vapply(S, am, 0, i = seq_len(nrow(d)))
  boot <- vapply(rep_idx, function(i) vapply(S, am, 0, i = i), numeric(ncol(S)))
  otros <- setdiff(names(S), "HAND")
  dif <- do.call(rbind, lapply(otros, function(k) {
    db <- boot["HAND", ] - boot[k, ]; db <- db[is.finite(db)]
    p <- min(1, 2 * min(mean(db <= 0), mean(db >= 0)))
    data.frame(contraste = paste("HAND -", k), dif_auc_media = round(obs[["HAND"]] - obs[[k]], 3),
               ic95_bloques = paste(round(stats::quantile(db, c(.025, .975)), 3), collapse = " to "),
               p_bilateral = p)
  }))
  dif$p_holm <- stats::p.adjust(dif$p_bilateral, "holm")
  dif$p_bilateral <- signif(dif$p_bilateral, 2); dif$p_holm <- signif(dif$p_holm, 2)
  mixtos <- do.call(rbind, lapply(names(S), function(k) {
    x <- data.frame(afectada = d$afectada, s = as.numeric(scale(S[[k]])), comuna = d$comuna)
    m <- suppressWarnings(lme4::glmer(afectada ~ s + (1 | comuna), x, family = stats::binomial()))
    data.frame(predictor = k, coef_estandarizado = round(lme4::fixef(m)[["s"]], 3), aic = round(stats::AIC(m), 1))
  }))
  mixtos$delta_aic <- mixtos$aic - min(mixtos$aic)
  list(auc_media = data.frame(predictor = names(obs), auc_media_6_comunas = round(obs, 3), n_gfi = n_gfi),
       diferencias = dif, mixtos = mixtos, bloques = length(unique(d$bloque)))
}

# 3. Barrido del peso de la amenaza en la agregación aditiva (unidades de Taltal).
barrido_aditivo <- function(d, pesos = seq(0, 1, 0.05)) {
  C <- (d$v_extendida_indice + d$exposicion + d$locc) / 3
  data.frame(peso_amenaza = pesos,
             auc_aditivo_alcance = vapply(pesos, function(w) round(auc(w * d$p_alcance + (1 - w) * C, d$afectada), 3), 0),
             auc_aditivo_susceptibilidad = vapply(pesos, function(w) round(auc(w * d$susceptibilidad + (1 - w) * C, d$afectada), 3), 0),
             auc_multiplicativo = round(auc(d$p_alcance * C, d$afectada), 3))
}
