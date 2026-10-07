# Evaluación del modelo de fricción simple de Flow-R (SFLM) propagado por D8 en las cuencas
# completas de las seis comunas (bin/sflm.py; fuentes por el criterio de eventos extremos de
# Horton et al. 2013). Para cada celda, "ángulo crítico" = mayor ángulo de alcance (grilla
# 0,25-15°) con el que la alcanza el flujo de alguna fuente; con tope de velocidad (15 m/s) y
# sin tope (ángulo de alcance H/L de la línea de energía). No alcanzada = 0.
#   - AUC por comuna (manzanas; mismas que la validación externa)
#   - ángulo requerido: mayor α con el que se alcanzaría el 50 / 80 % de las manzanas afectadas
#   - calibración en Atacama (LOCO, α que maximiza el índice de Youden con comunas de igual peso)
#     y transferencia a Taltal

region_comuna <- c("TALTAL" = "taltal", "CHAÑARAL" = "salado", "DIEGO DE ALMAGRO" = "salado",
                   "TIERRA AMARILLA" = "copiapo", "COPIAPO" = "copiapo", "ALTO DEL CARMEN" = "huasco")

bloques_sflm <- function(catastro, censo_gpkg, area_2015, dir_region, p, umbral = 2000) {
  do.call(rbind, lapply(seq_len(nrow(comunas_validacion)), function(i) {
    comuna <- comunas_validacion$comuna[i]
    m <- manzanas_con_distancia(comuna, comunas_validacion$cut[i], catastro, censo_gpkg, area_2015, p)
    if (comuna != "TALTAL") m <- m[m$dist_area <= 5000, ]
    r <- file.path(dir_region, region_comuna[[comuna]])
    dir_c <- dirname(salida(p, "validacion_externa", carpeta_ascii(comuna), "x"))
    x <- data.frame(comuna = comuna, afectada = m$afectada,
                    a15 = extraer(m, file.path(r, "angulo_v15.tif")),
                    ainf = extraer(m, file.path(r, "angulo_vinf.tif")),
                    hand = extraer(m, file.path(dir_c, sprintf("hand_costa_%d.tif", umbral))))
    x$a15[is.na(x$a15)] <- 0; x$ainf[is.na(x$ainf)] <- 0
    x
  }))
}

youden_alfa <- function(d, var, alfas = seq(0, 15, 0.25)) {
  w <- 1 / table(d$comuna)[d$comuna]
  j <- vapply(alfas, function(a) {
    alc <- d[[var]] >= a & d[[var]] > 0
    sens <- sum(w[alc & d$afectada]) / sum(w[d$afectada]); esp <- sum(w[!alc & !d$afectada]) / sum(w[!d$afectada])
    sens + esp - 1
  }, 0)
  alfas[which.max(j)]
}

evaluar_sflm <- function(d, semilla = 42) {
  cal <- setdiff(unique(d$comuna), "TALTAL")
  por_comuna <- do.call(rbind, lapply(split(d, d$comuna), function(x) {
    f <- function(v) { r <- auc_ic(v, x$afectada, semilla = semilla); round(c(r[["auc"]], r[["ic_inf"]], r[["ic_sup"]]), 3) }
    q <- function(v, pr) { a <- v[x$afectada]; if (!length(a)) NA else stats::quantile(a, 1 - pr, names = FALSE) }
    a15 <- f(x$a15); ai <- f(x$ainf); h <- f(-x$hand)
    data.frame(comuna = x$comuna[1], n = nrow(x), afectadas = sum(x$afectada),
               pct_afectadas_alcanzadas_11 = round(100 * mean(x$a15[x$afectada] >= 11), 1),
               pct_afectadas_alcanzadas_7 = round(100 * mean(x$a15[x$afectada] >= 7), 1),
               alfa50_v15 = q(x$a15, 0.5), alfa80_v15 = q(x$a15, 0.8),
               alfa50_inf = q(x$ainf, 0.5), alfa80_inf = q(x$ainf, 0.8),
               auc_v15 = a15[1], ic_v15 = paste(a15[2:3], collapse = "-"),
               auc_inf = ai[1], ic_inf = paste(ai[2:3], collapse = "-"),
               auc_hand = h[1], ic_hand = paste(h[2:3], collapse = "-"))
  }))
  loco <- do.call(rbind, lapply(c(cal, "TALTAL"), function(cm) {
    tr <- d[d$comuna %in% setdiff(cal, cm), ]; te <- d[d$comuna == cm, ]
    do.call(rbind, lapply(c("a15", "ainf"), function(v) {
      a <- youden_alfa(tr, v); alc <- te[[v]] >= a & te[[v]] > 0
      data.frame(comuna_prueba = cm, modelo = v, alfa_calibrado = a,
                 sensibilidad = round(mean(alc[te$afectada]), 3), especificidad = round(mean(!alc[!te$afectada]), 3))
    }))
  }))
  list(por_comuna = por_comuna, loco = loco)
}
