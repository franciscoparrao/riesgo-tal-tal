# Paquete de datos para el depósito en Zenodo (CC BY 4.0): tablas por unidad y manzana, índice de
# alcance a 30 m y tablas de resultados. No incluye la huella digitalizada de Taltal (derivada del
# mapa de Sernageomin; disponible a pedido), solo la clase afectada/no afectada por unidad.
exportar_deposito <- function(d, unidades, bloques, alcance_30m, csvs, p) {
  dir <- dirname(salida(p, "deposito", "x")); unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
  u <- sf::st_read(unidades, quiet = TRUE)[, "unidad"]
  d$riesgo_alcance <- d$p_alcance * (d$v_extendida_indice + d$exposicion + d$locc) / 3
  d$indice_aditivo_convencional <- 0.25 * (d$susceptibilidad + d$v_extendida_indice + d$exposicion + d$locc)
  ut <- merge(u, d[, c("unidad", "poblacion", "afectada", "frac_afectada", "susceptibilidad", "p_alcance",
                       "v_extendida_indice", "exposicion", "locc", "riesgo_alcance", "indice_aditivo_convencional")],
              by = "unidad")
  names(ut)[names(ut) == "p_alcance"] <- "indice_alcance"
  sf::st_write(ut, file.path(dir, "taltal_unidades_indices.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  b <- bloques[, c("comuna", "x", "y", "afectada", "hand", "dist", "elev", "pendiente", "area_cauce")]
  names(b) <- c("commune", "x_utm19s", "y_utm19s", "affected", "hand_m", "dist_channel_m", "elevation_m", "slope_deg", "receiving_channel_area_km2")
  utils::write.csv(b, file.path(dir, "blocks_six_communes.csv"), row.names = FALSE)
  terra::writeRaster(terra::rast(alcance_30m), file.path(dir, "reach_index_30m_taltal.tif"), overwrite = TRUE,
                     datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
  dir.create(file.path(dir, "results"))
  for (f in csvs) {   # etiqueta interna del índice aditivo -> nombre del paper
    x <- gsub("riesgo_aditivo_tesis", "indice_aditivo_convencional", readLines(f, encoding = "UTF-8"), fixed = TRUE)
    writeLines(x, file.path(dir, "results", basename(f)), useBytes = TRUE)
  }
  existe(dir)
}
