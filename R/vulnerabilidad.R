# Vulnerabilidad, modo corregido.
# Clustering no supervisado (K-means y HDBSCAN) sobre variables censales (Censo 2024) e índices
# espectrales urbanos (Sentinel-2 2024), con estas decisiones respecto de la versión convencional:
#  - sin x, y como variables; variables estandarizadas (z); semilla fija;
#  - k por coeficiente de silueta (Rousseeuw, 1987), como dice el texto;
#  - HDBSCAN: el ruido (cluster 0) NO es una clase: todo punto/celda se asigna al centroide
#    más cercano en el espacio estandarizado;
#  - reflectancia S2 corregida por BOA_ADD_OFFSET: (DN - 1000) / 10000 (issues_surtgis.md #9);
#  - VrNIRBI con su fórmula (no VgNIRBI).
# Niveles: regla automática (promedio de variables censales estandarizadas del centroide,
# orientadas "más = más vulnerable") con tabla de perfiles para revisión experta y posible
# sobrescritura en config.yml (vulnerabilidad$niveles_override).

# --- Variables censales por unidad -------------------------------------------------------

# Sin densidad poblacional: es exposición y dominaba el puntaje de los niveles
# (invertía urbano/rural según la variante).
variables_censales <- list(
  base = c("prop_0_13", "prop_60_mas", "sin_educacion", "sin_vivienda_propia", "migrantes"),
  extendida = c("prop_0_13", "prop_60_mas", "sin_educacion", "sin_vivienda_propia", "migrantes",
                "hacinamiento", "irrecuperables", "discapacidad", "agua_camion")
)

razon <- function(num, den) ifelse(is.na(den) | den <= 0, 0, pmin(num / den, 1))

unidades_vulnerabilidad <- function(censo, comuna_utm, epsg) {
  s0 <- function(...) Reduce(`+`, lapply(list(...), function(x) ifelse(is.na(x), 0, x)))
  cols <- c("n_per", "n_edad_0_5", "n_edad_6_13", "n_edad_60_mas", "n_inmigrantes",
            "n_cine_nunca_curso_primera_infancia", "n_cine_primaria", "n_cine_secundaria",
            "n_cine_terciaria_maestria_doctorado", "n_cine_especial_diferencial",
            "n_tenencia_propia_pagada", "n_tenencia_propia_pagandose",
            "n_tenencia_arrendada_contrato", "n_tenencia_arrendada_sin_contrato",
            "n_tenencia_cedida_trabajo", "n_tenencia_cedida_familiar", "n_tenencia_otro",
            "n_vp_ocupada", "n_viv_hacinadas", "n_viv_irrecuperables", "n_discapacidad",
            "n_fuente_agua_camion", "SHAPE")
  x <- do.call(rbind, lapply(censo$capas, function(capa) {
    sf::st_read(censo$ruta, quiet = TRUE, query = sprintf("SELECT %s FROM %s WHERE CUT=%d",
                                                          paste(cols, collapse = ", "), capa, censo$cut))
  }))
  x <- sf::st_transform(sf::st_make_valid(x), epsg)
  km2 <- as.numeric(sf::st_area(x)) / 1e6
  cine <- s0(x$n_cine_nunca_curso_primera_infancia, x$n_cine_primaria, x$n_cine_secundaria,
             x$n_cine_terciaria_maestria_doctorado, x$n_cine_especial_diferencial)
  tenencia_no_propia <- s0(x$n_tenencia_arrendada_contrato, x$n_tenencia_arrendada_sin_contrato,
                           x$n_tenencia_cedida_trabajo, x$n_tenencia_cedida_familiar,
                           x$n_tenencia_otro)
  tenencia <- tenencia_no_propia + s0(x$n_tenencia_propia_pagada, x$n_tenencia_propia_pagandose)
  u <- sf::st_sf(
    poblacion = s0(x$n_per),
    prop_0_13 = razon(s0(x$n_edad_0_5, x$n_edad_6_13), x$n_per),
    prop_60_mas = razon(s0(x$n_edad_60_mas), x$n_per),
    densidad_log = log1p(s0(x$n_per) / km2),
    sin_educacion = razon(s0(x$n_cine_nunca_curso_primera_infancia, x$n_cine_primaria), cine),
    sin_vivienda_propia = razon(tenencia_no_propia, tenencia),
    migrantes = razon(s0(x$n_inmigrantes), x$n_per),       # NA = secreto estadístico -> 0
    hacinamiento = razon(s0(x$n_viv_hacinadas), x$n_vp_ocupada),
    irrecuperables = razon(s0(x$n_viv_irrecuperables), x$n_vp_ocupada),
    discapacidad = razon(s0(x$n_discapacidad), x$n_per),
    agua_camion = razon(s0(x$n_fuente_agua_camion), x$n_vp_ocupada),
    geometry = sf::st_geometry(x)
  )
  # Resto de la comuna sin unidad censal: todas las variables en 0 (NaN -> 0).
  com <- sf::st_union(sf::st_read(comuna_utm, quiet = TRUE))
  u <- suppressWarnings(sf::st_intersection(u, com))
  resto <- sf::st_difference(com, sf::st_union(u))
  if (length(resto) > 0 && !all(sf::st_is_empty(resto))) {
    cero <- as.data.frame(lapply(sf::st_drop_geometry(u)[1, ], function(v) 0))
    u <- rbind(u, sf::st_sf(cero, geometry = resto))
  }
  u$unidad <- seq_len(nrow(u))
  u
}

unidades_archivo <- function(censo_ruta, comuna_utm, p) {
  destino <- salida(p, "componentes", "vulnerabilidad", "unidades_censales.gpkg")
  censo <- utils::modifyList(p$censo, list(ruta = censo_ruta))
  sf::st_write(unidades_vulnerabilidad(censo, comuna_utm, p$crs), destino,
               delete_dsn = TRUE, quiet = TRUE)
  existe(destino)
}

# --- Índices espectrales -----------------------------------------------------------------

indices_nombres <- c("blfei", "brba", "nbai", "nhfd", "ui", "vgnirbi", "vrnirbi", "pisi")

# b: data.frame/lista con B01..B12 en DN del composite; devuelve los 8 índices.
indices_espectrales <- function(b) {
  r <- lapply(b, function(x) (x - 1000) / 10000)          # BOA_ADD_OFFSET (issue #9)
  nd <- function(a, c) (a - c) / (a + c)
  m <- (r$B03 + r$B04 + r$B12) / 3
  data.frame(
    blfei = nd(m, r$B11),
    brba = r$B04 / r$B11,
    nbai = nd(r$B12, r$B11),                                # ((B12-B11)/B3)/((B12+B11)/B3)
    nhfd = nd(r$B05, r$B01),
    ui = nd(r$B12, r$B08),
    vgnirbi = nd(r$B03, r$B08),
    vrnirbi = nd(r$B04, r$B08),
    pisi = 0.8192 * r$B02 - 0.5735 * r$B08 + 0.0750
  )
}

bandas_s2 <- c("B01", "B02", "B03", "B04", "B05", "B08", "B11", "B12")

# --- Muestra de entrenamiento ------------------------------------------------------------

muestra_vulnerabilidad <- function(unidades, dir_s2, p) {
  u <- sf::st_read(unidades, quiet = TRUE)
  set.seed(p$semilla)
  pts <- sf::st_sample(sf::st_union(u), p$vulnerabilidad$n_muestra)
  pts <- sf::st_sf(geometry = pts)
  censo <- sf::st_drop_geometry(sf::st_join(pts, u, join = sf::st_intersects, left = TRUE))
  censo <- censo[!duplicated(seq_len(nrow(pts))), , drop = FALSE]
  b <- terra::rast(file.path(dir_s2, paste0(bandas_s2, ".tif")))
  names(b) <- bandas_s2
  bv <- terra::extract(b, terra::vect(pts), ID = FALSE)
  xy <- sf::st_coordinates(pts)
  d <- cbind(x = xy[, 1], y = xy[, 2], censo, indices_espectrales(bv))
  d <- d[stats::complete.cases(d) & is.finite(rowSums(d[, indices_nombres])), ]
  d
}

# --- Espectral por unidad ---------------------------------------------------------------

# Muestra regular de píxeles dentro de cada unidad poblada (1 cada `paso_m2`, entre 1 y
# `max_px`), 8 índices por píxel, PCA sobre todos los píxeles (estandarizado) y promedio de los
# primeros `n_pc` componentes por unidad. Reduce la redundancia del bloque espectral (PC1-3
# explican ~80 %) y evita que pese 8 variables de 18.
espectral_unidades <- function(unidades, dir_s2, p) {
  pv <- p$vulnerabilidad
  u <- sf::st_read(unidades, quiet = TRUE)
  u <- u[u$poblacion > 0, ]
  n_px <- pmin(pmax(ceiling(as.numeric(sf::st_area(u)) / pv$espectral_paso_m2), 1), pv$espectral_max_px)
  set.seed(p$semilla)
  pts <- lapply(seq_len(nrow(u)), function(i) {
    g <- sf::st_sample(sf::st_geometry(u)[i], n_px[i], type = "regular")
    if (length(g) == 0) g <- sf::st_point_on_surface(sf::st_geometry(u)[i])
    sf::st_sf(unidad = u$unidad[i], geometry = g)
  })
  pts <- do.call(rbind, pts)
  b <- terra::rast(file.path(dir_s2, paste0(bandas_s2, ".tif")))
  names(b) <- bandas_s2
  px <- cbind(unidad = pts$unidad,
              indices_espectrales(terra::extract(b, terra::vect(pts), ID = FALSE)))
  px <- px[stats::complete.cases(px) & is.finite(rowSums(px[indices_nombres])), ]
  pca <- stats::prcomp(px[indices_nombres], scale. = TRUE)
  pcs <- paste0("pc", seq_len(pv$n_pc))
  sc <- as.data.frame(pca$x[, seq_len(pv$n_pc), drop = FALSE]); names(sc) <- pcs
  tabla <- stats::aggregate(sc, list(unidad = px$unidad), mean)
  tabla$n_pixeles <- as.vector(table(px$unidad)[as.character(tabla$unidad)])
  list(tabla = tabla, pcs = pcs,
       varianza = pca$sdev^2 / sum(pca$sdev^2),
       cargas = round(pca$rotation[, seq_len(pv$n_pc), drop = FALSE], 3))
}

# --- Clustering de unidades --------------------------------------------------------------

estandarizar <- function(d, vars) {
  mu <- vapply(d[vars], mean, 0)
  sdv <- vapply(d[vars], stats::sd, 0)
  sdv[sdv == 0] <- 1
  list(z = scale(d[vars], center = mu, scale = sdv), mu = mu, sd = sdv)
}

silueta_k <- function(z, ks, semilla) {
  dz <- stats::dist(z)
  vapply(ks, function(k) {
    set.seed(semilla)
    km <- stats::kmeans(z, centers = k, nstart = 50, iter.max = 100)
    mean(cluster::silhouette(km$cluster, dz)[, "sil_width"])
  }, 0)
}

centroides <- function(z, grupo) {
  g <- sort(unique(grupo))
  t(vapply(g, function(k) colMeans(z[grupo == k, , drop = FALSE]), numeric(ncol(z))))
}

asignar_centroide <- function(z, cent) {
  d2 <- vapply(seq_len(nrow(cent)), function(k) colSums((t(z) - cent[k, ])^2), numeric(nrow(z)))
  if (is.null(dim(d2))) d2 <- matrix(d2, nrow = 1)
  max.col(-d2, ties.method = "first")
}

# Niveles: promedio de variables censales estandarizadas del centroide (todas orientadas
# "más = más vulnerable"); rango -> 0.15, ..., 1.00 (7 niveles si k = 7).
niveles_automaticos <- function(cent, vars_censo) {
  puntaje <- rowMeans(cent[, vars_censo, drop = FALSE])
  k <- nrow(cent)
  niveles <- seq(0.15, 1.00, length.out = k)
  data.frame(cluster = seq_len(k), puntaje = puntaje, nivel = niveles[rank(puntaje, ties.method = "first")])
}

# Fuzzificador m recomendado según N observaciones y D dimensiones (Schwämmle & Jensen, 2010,
# Bioinformatics 26:2841). Con m = 2 y D = 8-12 las pertenencias tienden a ser uniformes y la
# vulnerabilidad degenera a una constante (observado: 0.57-0.58 en la variante extendida).
m_fuzzy <- function(N, D) {
  1 + (1418 / N + 22.05) * D^-2 + (12.33 / N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
}

# Unidad de análisis = unidad censal poblada (manzana urbana o entidad rural), una observación
# por unidad; la población por cluster se reporta en el perfil. Fuera de ellas: no aplica.
modelo_vulnerabilidad <- function(unidades, espectral, variante, metodo, p) {
  pv <- p$vulnerabilidad
  u <- sf::st_drop_geometry(sf::st_read(unidades, quiet = TRUE))
  d <- merge(u[u$poblacion > 0, ], espectral$tabla, by = "unidad", all.x = TRUE)
  d[espectral$pcs][is.na(d[espectral$pcs])] <- 0
  vars_censo <- variables_censales[[variante]]
  vars <- c(vars_censo, espectral$pcs)
  e <- estandarizar(d, vars)
  z <- e$z
  ks <- pv$k_rango[1]:min(pv$k_rango[2], nrow(z) - 1)
  sil <- silueta_k(z, ks, p$semilla)
  k <- if (identical(pv$k, "auto")) ks[which.max(sil)] else as.integer(pv$k)
  set.seed(p$semilla)
  info <- list(minPts = NA, ruido = 0, m = NA)
  if (metodo == "indice") {
    # Índice continuo (estilo SoVI): promedio de variables censales estandarizadas, 0-1.
    puntaje <- rowMeans(z[, vars_censo, drop = FALSE])
    d$nivel <- (puntaje - min(puntaje)) / (max(puntaje) - min(puntaje))
    d$cluster <- as.integer(cut(d$nivel, stats::quantile(d$nivel, 0:5 / 5), include.lowest = TRUE))
    cent <- centroides(z, d$cluster)          # quintiles del índice, solo para el perfil
    niveles <- data.frame(cluster = 1:5, puntaje = NA_real_,
                          nivel = as.vector(tapply(d$nivel, d$cluster, mean)))
  } else if (metodo == "fuzzy") {
    m_f <- if (identical(pv$fuzzy_m, "auto")) m_fuzzy(nrow(z), ncol(z)) else as.numeric(pv$fuzzy_m)
    fc <- e1071::cmeans(z, centers = k, m = m_f, iter.max = 500, method = "cmeans")
    cent <- fc$centers
    niveles <- niveles_automaticos(cent, vars_censo)
    override <- pv$niveles_override[[paste(variante, metodo, sep = "_")]]
    if (!is.null(override)) niveles$nivel[as.integer(names(override))] <- unlist(override)
    d$nivel <- as.vector(fc$membership %*% niveles$nivel)   # gradiente: sum(u_ik * nivel_k)
    d$cluster <- fc$cluster
    info$m <- m_f
  } else {
    if (metodo == "kmeans") {
      cent <- stats::kmeans(z, centers = k, nstart = 50, iter.max = 100)$centers
    } else {
      mejor <- NULL
      for (mp in pv$hdbscan_minpts) {
        h <- dbscan::hdbscan(z, minPts = mp)
        nk <- length(setdiff(unique(h$cluster), 0))
        if (nk >= 2 && (is.null(mejor) || abs(nk - k) < abs(mejor$nk - k))) {
          mejor <- list(h = h, nk = nk, mp = mp)
        }
        if (!is.null(mejor) && mejor$nk == k) break
      }
      if (is.null(mejor)) stop("HDBSCAN no encontró >= 2 clusters con los minPts dados")
      cl <- mejor$h$cluster
      cent <- centroides(z[cl > 0, , drop = FALSE], cl[cl > 0])
      info$minPts <- mejor$mp; info$ruido <- mean(cl == 0)
    }
    d$cluster <- asignar_centroide(z, cent)        # ruido de HDBSCAN -> centroide más cercano
    niveles <- niveles_automaticos(cent, vars_censo)
    override <- pv$niveles_override[[paste(variante, metodo, sep = "_")]]
    if (!is.null(override)) niveles$nivel[as.integer(names(override))] <- unlist(override)
    d$nivel <- niveles$nivel[d$cluster]
  }
  perfil <- stats::aggregate(d[vars], list(cluster = d$cluster), mean)
  tam <- stats::aggregate(cbind(n_unidades = 1, poblacion = d$poblacion), list(cluster = d$cluster), sum)
  tam$pct_poblacion <- 100 * tam$poblacion / sum(tam$poblacion)
  list(variante = variante, metodo = metodo, vars = vars, vars_censo = vars_censo,
       mu = e$mu, sd = e$sd, centroides = cent, k = nrow(cent), k_silueta = k,
       silueta = data.frame(k = ks, silueta = sil), info = info,
       asignacion = d[, c("unidad", "cluster", "nivel")],
       niveles = merge(merge(niveles, tam, by = "cluster"), perfil, by = "cluster"))
}

# Tabla de perfiles para revisión experta (como Tablas 4.3/4.4).
perfiles_csv <- function(modelo, p) {
  destino <- salida(p, "componentes", "vulnerabilidad",
                    sprintf("perfiles_%s_%s.csv", modelo$variante, modelo$metodo))
  n <- modelo$niveles
  n[] <- lapply(n, function(x) if (is.numeric(x)) round(x, 4) else x)
  utils::write.csv(n[order(-n$nivel), ], destino, row.names = FALSE)
  existe(destino)
}

# --- Raster ------------------------------------------------------------------------------

# Nivel por unidad poblada -> 10 m con SurtGIS. Fuera de unidades pobladas: NaN (no aplica).
vulnerabilidad_raster <- function(modelo, unidades, plantilla_fina, p) {
  nombre <- paste(modelo$variante, modelo$metodo, sep = "_")
  vec <- salida(p, "componentes", "vulnerabilidad", paste0("vulnerabilidad_", nombre, ".gpkg"))
  destino <- salida(p, "componentes", paste0("vulnerabilidad_", nombre, ".tif"))
  u <- sf::st_read(unidades, quiet = TRUE)
  u <- merge(u[, "unidad"], modelo$asignacion, by = "unidad")
  sf::st_write(u, vec, delete_dsn = TRUE, quiet = TRUE)
  surtgis("rasterize", "--reference", plantilla_fina, "--attribute", "nivel", vec, destino)
  existe(destino)
}

# Relleno de huecos (camanchaca costera enmascarada por SCL en las 4 fechas del composite
# principal): cada celda vacía toma el valor del composite costero (todo 2024, hasta 12
# escenas), alineado a la misma grilla de 20 m. Las celdas con dato no se tocan.
rellenar_s2 <- function(bandas, bandas_relleno, p) {
  terra::terraOptions(memfrac = 0.2, progress = 0)
  vapply(seq_along(bandas), function(i) {
    destino <- salida(p, "s2_2024", basename(bandas[i]))
    r <- terra::rast(bandas[i])
    f <- terra::extend(terra::rast(bandas_relleno[i]), r)
    terra::cover(r, f, filename = destino, overwrite = TRUE,
                 wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")))
    existe(destino)
  }, "")
}
