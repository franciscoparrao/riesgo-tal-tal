# LOCC binaria a partir de una isócrona (usada por R/locc.R).
locc_poligonos <- function(comuna, isocrona, cfg, archivo = "locc_poligonos.gpkg") {
  destino <- salida(cfg, "tmp", archivo)
  iso <- sf::st_union(sf::st_zm(sf::st_read(isocrona, quiet = TRUE)))
  com <- sf::st_union(sf::st_transform(sf::st_read(comuna, quiet = TRUE), sf::st_crs(iso)))
  p <- sf::st_sf(locc = c(0, 1),
                 geometry = c(sf::st_intersection(com, iso), sf::st_difference(com, iso)))
  sf::st_write(p, destino, delete_dsn = TRUE, quiet = TRUE)
  existe(destino)
}
