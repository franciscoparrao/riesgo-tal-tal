# Ángulo de alcance requerido para que un flujo llegue a Taltal (diagnóstico de Flow-Py/Flow-R).
# Con un modelo de fricción simple (línea de energía de pendiente tan α) y tope de velocidad
# (altura de energía <= v_max^2 / 2g, 11,5 m para 15 m/s), un flujo se detiene a pocas decenas de
# metros de donde la pendiente del cauce cae bajo α. Se calcula, con D8:
#   (1) para cada celda fuente (susceptibilidad >= cuantil q) de la cuenca de la Quebrada de
#       Taltal, el ángulo de la línea fuente -> entrada del área afectada al mapa, sin tope de
#       velocidad: atan((z_fuente - z_entrada) / L_recorrido);
#   (2) la pendiente media del cauce entre la entrada al mapa y el mar y por tramos de 1 km.

angulo_requerido <- function(amenaza, marco, p, q = 0.99, entrada = c(353917.3, 7185196)) {
  stopifnot(length(amenaza) > 0)
  dir_a <- dirname(salida(p, "amenaza", "x"))
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = dirname(salida(p, "tmp", "terra", "x")))
  fd <- terra::rast(file.path(dir_a, "flowdir_terra.tif"))
  dem <- terra::rast(file.path(dir_a, "dem_fill.tif"))
  s <- terra::rast(file.path(dir_a, "susc_30m.tif"))
  ws <- terra::classify(terra::watershed(fd, rbind(entrada)), cbind(0, NA))
  e <- terra::ext(terra::trim(ws)) + 5000
  e <- terra::union(e, terra::ext(sf::st_read(marco, quiet = TRUE)) + 2000)
  fd <- terra::crop(fd, e); dem <- terra::crop(dem, e); s <- terra::crop(s, e); ws <- terra::crop(ws, e)
  nc <- terra::ncol(fd); nr <- terra::nrow(fd)
  dv <- terra::values(fd, mat = FALSE); z <- terra::values(dem, mat = FALSE)
  # terra flowdir: 1=E 2=SE 4=S 8=SW 16=W 32=NW 64=N 128=NE
  cod <- c(1, 2, 4, 8, 16, 32, 64, 128)
  dc <- c(1, 1, 0, -1, -1, -1, 0, 1); dr <- c(0, 1, 1, 1, 0, -1, -1, -1)
  paso <- 30 * c(1, sqrt(2), 1, sqrt(2), 1, sqrt(2), 1, sqrt(2))
  k <- match(dv, cod)
  i <- seq_along(dv) - 1L; r0 <- i %/% nc; c0 <- i %% nc
  r1 <- r0 + dr[k]; c1 <- c0 + dc[k]
  aguas_abajo <- ifelse(!is.na(k) & r1 >= 0 & r1 < nr & c1 >= 0 & c1 < nc, r1 * nc + c1 + 1L, NA_integer_)
  largo <- ifelse(is.na(k), NA_real_, paso[k])
  # Distancia a lo largo del flujo hasta la entrada, solo dentro de la cuenca, procesando de aguas
  # abajo a aguas arriba (acumulación decreciente).
  ce <- terra::cellFromXY(fd, rbind(entrada))
  en_cuenca <- which(!is.na(terra::values(ws, mat = FALSE)))
  acc <- terra::values(terra::flowAccumulation(fd), mat = FALSE)
  orden <- en_cuenca[order(acc[en_cuenca], decreasing = TRUE)]
  dist <- rep(NA_real_, length(dv)); dist[ce] <- 0
  for (j in orden) if (j != ce) { d <- aguas_abajo[j]; if (!is.na(d) && !is.na(dist[d])) dist[j] <- dist[d] + largo[j] }
  sv <- terra::values(s, mat = FALSE)
  umbral <- stats::quantile(sv[en_cuenca], q, na.rm = TRUE)
  fuentes <- en_cuenca[!is.na(sv[en_cuenca]) & sv[en_cuenca] >= umbral]
  L <- dist[fuentes]
  ang <- atan((z[fuentes] - z[ce]) / L) * 180 / pi
  # Perfil del cauce desde la entrada hasta el mar
  perfil <- c(); lp <- c(); j <- ce; l <- 0
  while (!is.na(j) && length(perfil) < 1e5) { perfil <- c(perfil, j); lp <- c(lp, l); l <- l + largo[j]; j <- aguas_abajo[j] }
  zp <- z[perfil]
  mar <- which(zp <= 1)[1]; if (!is.na(mar)) { zp <- zp[seq_len(mar)]; lp <- lp[seq_len(mar)] }
  tramos <- cut(lp, seq(0, max(lp) + 1000, 1000), include.lowest = TRUE)
  pend <- tapply(seq_along(zp), tramos, function(ii) atan((zp[min(ii)] - zp[max(ii)]) /
                                                      max(lp[max(ii)] - lp[min(ii)], 1)) * 180 / pi)
  resumen <- data.frame(
    metrica = c("n fuentes (cuenca, S >= q99)", "angulo requerido fuente->entrada, mediana (grados)",
                "angulo requerido, p90", "angulo requerido, maximo",
                "% fuentes que requieren < 7 grados", "% fuentes que requieren < 11 grados",
                "cota entrada (m)", "largo entrada->mar (m)", "pendiente media entrada->mar (grados)",
                "pendiente maxima por tramo de 1 km (grados)"),
    valor = round(c(length(fuentes), stats::median(ang, na.rm = TRUE), stats::quantile(ang, 0.9, na.rm = TRUE),
                    max(ang, na.rm = TRUE), 100 * mean(ang < 7, na.rm = TRUE), 100 * mean(ang < 11, na.rm = TRUE),
                    z[ce], max(lp), atan(zp[1] / max(lp)) * 180 / pi, max(pend, na.rm = TRUE)), 2))
  list(resumen = resumen, angulos = data.frame(angulo = round(ang[is.finite(ang)], 3)))
}
