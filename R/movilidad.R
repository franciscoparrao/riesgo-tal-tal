# Movilidad de los aluviones de marzo de 2015: pendiente de los cauces por los que pasó y se
# depositó el flujo en las seis comunas (áreas afectadas de Sernageomin), para compararla con los
# ángulos mínimos de alcance publicados para flujos de detritos (Horton et al. 2013: 11° gruesos,
# 7° finos). Mismo DEM (Copernicus GLO-30, 30 m, rellenado) y red de cauces (>= umbral celdas) que
# la validación externa.
#   pendiente de cauce = atan(Δz / L) entre cada celda de cauce afectada y el punto L = `largo` m
#   aguas abajo por D8 (o hasta el mar / borde si está más cerca, con L >= 150 m)

pendiente_cauce_afectado <- function(dir, a, umbral = 2000, largo = 500, margen = 3000) {
  dem <- terra::rast(file.path(dir, "dem_fill.tif"))
  e <- terra::intersect(terra::ext(terra::vect(a)) + margen, terra::ext(dem))
  dem <- terra::crop(dem, e)
  acc <- terra::crop(terra::rast(file.path(dir, "facc.tif")), e)
  fd <- terra::terrain(dem, v = "flowdir")
  afect <- terra::rasterize(terra::vect(a), dem, touches = TRUE)
  cauce <- which(terra::values(acc, mat = FALSE) > umbral & !is.na(terra::values(afect, mat = FALSE)))
  z <- terra::values(dem, mat = FALSE); dv <- terra::values(fd, mat = FALSE)
  nc <- terra::ncol(dem); nr <- terra::nrow(dem)
  cod <- c(1, 2, 4, 8, 16, 32, 64, 128)
  dc <- c(1, 1, 0, -1, -1, -1, 0, 1); dr <- c(0, 1, 1, 1, 0, -1, -1, -1)
  paso <- 30 * c(1, sqrt(2), 1, sqrt(2), 1, sqrt(2), 1, sqrt(2))
  pos <- cauce; L <- numeric(length(pos)); activo <- rep(TRUE, length(pos))
  for (it in seq_len(ceiling(largo / 30) + 5)) {
    k <- match(dv[pos], cod)
    i <- pos - 1L; r1 <- i %/% nc + dr[k]; c1 <- i %% nc + dc[k]
    ok <- activo & !is.na(k) & r1 >= 0 & r1 < nr & c1 >= 0 & c1 < nc
    nuevo <- ifelse(ok, r1 * nc + c1 + 1L, pos)
    ok <- ok & !is.na(z[nuevo])
    L[ok] <- L[ok] + paso[k[ok]]; pos[ok] <- nuevo[ok]
    activo <- ok & L < largo
    if (!any(activo)) break
  }
  g <- ifelse(L >= 150, atan(pmax(z[cauce] - z[pos], 0) / L) * 180 / pi, NA_real_)
  pend <- terra::values(terra::terrain(dem, "slope", unit = "degrees"), mat = FALSE)
  list(gradiente = g[is.finite(g)], pendiente_area = pend[!is.na(terra::values(afect, mat = FALSE))])
}

movilidad_comunas <- function(catastro, area_2015, manzanas_validacion, p, umbral = 2000, largo = 500) {
  stopifnot(nrow(manzanas_validacion) > 0)     # dependencia: DEM por comuna ya descargado
  terra::terraOptions(memfrac = 0.1, progress = 0, tempdir = dirname(salida(p, "tmp", "terra", "x")))
  res <- lapply(comunas_validacion$comuna, function(comuna) {
    a <- sf::st_sf(geometry = afectada_comuna(catastro, comuna, p$crs, area_2015))
    dir <- dirname(salida(p, "validacion_externa", carpeta_ascii(comuna), "x"))
    x <- pendiente_cauce_afectado(dir, a, umbral, largo)
    list(comuna = comuna, g = x$gradiente, s = x$pendiente_area)
  })
  resumen <- do.call(rbind, lapply(res, function(r) {
    q <- function(v, pr) unname(round(stats::quantile(v, pr, na.rm = TRUE), 2))
    data.frame(comuna = r$comuna, celdas_cauce = length(r$g), largo_cauce_km = round(length(r$g) * 0.03, 1),
               gradiente_p10 = q(r$g, 0.1), gradiente_mediana = q(r$g, 0.5), gradiente_p90 = q(r$g, 0.9),
               pct_bajo_3 = round(100 * mean(r$g < 3), 1), pct_bajo_7 = round(100 * mean(r$g < 7), 1),
               pct_bajo_11 = round(100 * mean(r$g < 11), 1), pendiente_area_mediana = q(r$s, 0.5))
  }))
  celdas <- do.call(rbind, lapply(res, function(r) data.frame(comuna = r$comuna, gradiente = r$g)))
  list(resumen = resumen, celdas = celdas)
}
