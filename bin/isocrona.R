#!/usr/bin/env Rscript
# Isócrona de un tipo de infraestructura, en un proceso aislado (lo lanza isocrona_tipo()).
# Uso: Rscript bin/isocrona.R <network_dir> <instalaciones.gpkg> <umbral_min> <destino.gpkg> <xmx> <lote>
# Procesa los orígenes en lotes y va uniendo los polígonos: la memoria no crece con el número
# de orígenes (526 puntos para las líneas de transmisión muestreadas cada 500 m).

a <- commandArgs(trailingOnly = TRUE)
red_dir <- a[1]; instalaciones <- a[2]; umbral <- as.numeric(a[3]); destino <- a[4]
options(java.parameters = paste0("-Xmx", a[5]))
lote <- as.integer(a[6])

suppressPackageStartupMessages({ library(sf); library(r5r) })
p <- st_read(instalaciones, quiet = TRUE)
net <- build_network(red_dir, verbose = FALSE)

acumulado <- NULL
grupos <- split(seq_len(nrow(p)), ceiling(seq_len(nrow(p)) / lote))
for (k in seq_along(grupos)) {
  g <- grupos[[k]]
  iso <- isochrone(net, origins = p[g, ], mode = "CAR", cutoffs = umbral,
                   polygon_output = TRUE, sample_size = 1,
                   departure_datetime = as.POSIXct("2026-09-29 10:00:00"),
                   max_trip_duration = umbral + 5, progress = FALSE)
  iso <- iso[iso$isochrone == umbral & !st_is_empty(iso), ]
  if (nrow(iso) > 0) {
    u <- st_union(st_make_valid(st_transform(iso, 32719)))
    acumulado <- if (is.null(acumulado)) u else st_union(acumulado, u)
  }
  rm(iso); invisible(gc())
  cat(sprintf("lote %d/%d  RSS R=%.0f MB\n", k, length(grupos),
              as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024))
}
stop_r5(net)

if (is.null(acumulado)) stop("Ningún origen alcanzó la red vial")
st_write(st_sf(umbral_min = umbral, n_origenes = nrow(p), geometry = acumulado),
         destino, delete_dsn = TRUE, quiet = TRUE)
cat("OK", destino, "\n")
