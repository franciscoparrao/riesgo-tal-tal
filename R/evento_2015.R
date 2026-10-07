# Evento de validación: aluvión del 25-26 de marzo de 2015 en Taltal.
# Fuente: Sernageomin (con ONEMI, CIREN, SAF, DigitalGlobe), "Zonas fotointerpretadas
# (controladas en terreno)", fotointerpretación sobre imagen SAF 2015-03-27, Google Earth y
# cartografía geológica de Sernageomin; escala 1:20.000; versión del 23 de abril de 2015.
# Solo se dispone de dos capturas de pantalla del mapa (norte y sur). Se georreferencian con la
# grilla UTM de 1 km (líneas detectadas automáticamente, ~4,6 m/px) y se digitaliza por color
# la clase "Área afectada" (roja). Precisión esperable: ±10-20 m (declarar en el paper).

# Líneas de grilla detectadas (píxel, 1-based) -> coordenadas UTM (EPSG:32719, a verificar).
grilla_2015 <- list(
  norte = list(x_px = c(260, 476, 693, 910, 1126), x_utm = seq(350000, 354000, 1000),
               y_px = c(27, 243, 460, 677), y_utm = c(7191000, 7190000, 7189000, 7188000)),
  sur = list(x_px = c(260, 476, 693, 910, 1126), x_utm = seq(350000, 354000, 1000),
             y_px = c(115, 333, 549, 765), y_utm = c(7187000, 7186000, 7185000, 7184000))
)
# Elementos cartográficos rojos que no son datos (xmin, xmax, ymin, ymax en píxeles):
# logo "Gobierno de Chile" (norte); leyenda con muestra "Área afectada" y logo ONEMI (sur).
excluir_2015 <- list(
  norte = list(c(55, 245, 0, 155)),
  sur = list(c(55, 350, 440, 840), c(345, 505, 455, 790))
)

# Raster RGB georreferenciado (ajuste lineal por eje; centro del píxel i = i - 0.5).
georref_captura <- function(png_ruta, g, epsg, excluir = list()) {
  m <- png::readPNG(png_ruta)[, , 1:3] * 255
  for (b in excluir) m[b[3]:b[4], b[1]:b[2], ] <- 0     # negro: no pasa el umbral rojo
  fx <- stats::lm(g$x_utm ~ I(g$x_px - 0.5)); fy <- stats::lm(g$y_utm ~ I(g$y_px - 0.5))
  bx <- stats::coef(fx); by <- stats::coef(fy)
  r <- terra::rast(m)
  terra::ext(r) <- c(bx[1], bx[1] + bx[2] * ncol(r), by[1] + by[2] * nrow(r), by[1])
  terra::crs(r) <- paste0("EPSG:", epsg)
  names(r) <- c("R", "G", "B")
  attr(r, "residuo_m") <- max(abs(c(stats::residuals(fx), stats::residuals(fy))))
  r
}

# Máscara "Área afectada": rojo saturado y claro (relleno rosado semitransparente y contorno
# rojo). Excluye "acopio" (granate oscuro), "evacuación" (amarillo) y "campamento" (verde).
mascara_afectada <- function(r) {
  R <- r$R; G <- r$G; B <- r$B
  roja <- (R > 170) & (R - G > 55) & (R - B > 55)
  # cierre morfológico para unir el achurado y los contornos de los tramos
  roja <- terra::focal(roja, w = 5, fun = "max", na.rm = TRUE)
  roja <- terra::focal(roja, w = 5, fun = "min", na.rm = TRUE)
  terra::classify(roja, cbind(0, NA))
}

digitalizar_2015 <- function(png_norte, png_sur, p, area_min_m2 = 2000) {
  destino <- salida(p, "validacion", "area_afectada_2015.gpkg")
  partes <- lapply(list(norte = png_norte, sur = png_sur), function(f) {
    nombre <- if (identical(f, png_norte)) "norte" else "sur"
    r <- georref_captura(f, grilla_2015[[nombre]], p$crs, excluir_2015[[nombre]])
    v <- terra::as.polygons(mascara_afectada(r), dissolve = TRUE)
    sf::st_as_sf(v)
  })
  a <- sf::st_union(do.call(rbind, partes))
  a <- sf::st_cast(sf::st_make_valid(a), "POLYGON")
  a <- sf::st_sf(geometry = a)
  a <- a[as.numeric(sf::st_area(a)) >= area_min_m2, ]
  # rellenar huecos internos (achurado) y suavizar el borde de píxel
  a$geometry <- sf::st_geometry(a) |> lapply(function(g) sf::st_polygon(g[1])) |> sf::st_sfc(crs = sf::st_crs(a))
  a$area_m2 <- as.numeric(sf::st_area(a))
  a$fuente <- "Sernageomin 2015 (v. 23-04-2015), digitalizado de captura"
  sf::st_write(a, destino, delete_dsn = TRUE, quiet = TRUE)
  existe(destino)
}
