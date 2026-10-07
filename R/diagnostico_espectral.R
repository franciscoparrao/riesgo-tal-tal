# Diagnóstico: ¿los índices espectrales urbanos capturan vulnerabilidad o geología?
# Sobre la muestra de vulnerabilidad (puntos aleatorios en la comuna):
#  1. Correlación (Spearman) entre los 8 índices y varianza explicada por PCA del bloque.
#  2. Tres clusterings K-means con el mismo k: censo + espectral (índice convencional), solo censo y
#     solo espectral.
#  3. Acuerdo (NMI, ARI) de cada clustering con las unidades geológicas (Sernageomin 1:1M) y
#     entre clusterings; línea base = etiquetas permutadas. En toda la comuna y solo en puntos
#     poblados (densidad > 0).

diagnostico_espectral <- function(muestra, geologia, p, k = 7, variante = "extendida") {
  dir <- salida(p, "diagnostico", "x")
  dir <- dirname(dir)
  d <- muestra
  pts <- sf::st_as_sf(d[, c("x", "y")], coords = c("x", "y"), crs = p$crs)
  geo <- sf::st_transform(sf::st_make_valid(sf::st_read(geologia, quiet = TRUE)), p$crs)
  idx <- sf::st_intersects(pts, geo)
  d$geologia <- geo$cod_geologia[vapply(idx, function(i) if (length(i)) i[1] else NA_integer_, 1L)]
  cobertura <- mean(!is.na(d$geologia))   # la capa recortada no cubre toda la comuna
  d <- d[!is.na(d$geologia), ]

  # 1. Correlación y PCA del bloque espectral
  corr <- stats::cor(d[indices_nombres], method = "spearman")
  utils::write.csv(round(corr, 3), file.path(dir, "correlacion_indices.csv"))
  pca <- stats::prcomp(d[indices_nombres], scale. = TRUE)
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  utils::write.csv(data.frame(componente = paste0("PC", seq_along(var_exp)),
                              varianza = round(var_exp, 4), acumulada = round(cumsum(var_exp), 4)),
                   file.path(dir, "pca_indices.csv"), row.names = FALSE)

  # 2. Clusterings
  vc <- variables_censales[[variante]]
  conjuntos <- list(censo_espectral = c(vc, indices_nombres), solo_censo = vc,
                    solo_espectral = indices_nombres)
  etiquetas <- lapply(conjuntos, function(v) {
    z <- scale(d[v]); z[!is.finite(z)] <- 0
    set.seed(p$semilla)
    stats::kmeans(z, centers = k, nstart = 25, iter.max = 200, algorithm = "MacQueen")$cluster
  })
  set.seed(p$semilla)
  etiquetas$aleatorio <- sample(etiquetas$censo_espectral)

  # 3. Acuerdo
  acuerdo <- function(a, b) {
    a <- as.integer(factor(a)); b <- as.integer(factor(b))   # aricode no acepta texto
    c(NMI = aricode::NMI(a, b), ARI = aricode::ARI(a, b))
  }
  filas <- list()
  for (ambito in c("comuna", "poblado")) {
    s <- if (ambito == "comuna") rep(TRUE, nrow(d)) else d$densidad_log > 0
    for (n in names(etiquetas)) {
      g <- acuerdo(etiquetas[[n]][s], d$geologia[s])
      cs <- acuerdo(etiquetas[[n]][s], etiquetas$solo_censo[s])
      filas[[length(filas) + 1]] <- data.frame(
        ambito = ambito, n_puntos = sum(s), clustering = n,
        NMI_geologia = g[["NMI"]], ARI_geologia = g[["ARI"]],
        NMI_con_solo_censo = cs[["NMI"]], ARI_con_solo_censo = cs[["ARI"]])
    }
  }
  tabla <- do.call(rbind, filas)
  tabla$cobertura_geologia <- round(cobertura, 3)
  tabla[, -(1:3)] <- round(tabla[, -(1:3)], 4)
  destino <- file.path(dir, "acuerdo_clusterings.csv")
  utils::write.csv(tabla, destino, row.names = FALSE)
  existe(destino)
}
