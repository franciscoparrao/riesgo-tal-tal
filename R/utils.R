# Utilidades: configuración, rutas de salida y llamada a SurtGIS.

leer_config <- function(modo = Sys.getenv("RIESGO_MODO", "corregido"), archivo = "config.yml") {
  todo <- yaml::read_yaml(archivo)
  if (!modo %in% names(todo)) stop("Modo desconocido: ", modo)
  cfg <- utils::modifyList(todo$default, todo[[modo]])
  cfg$modo <- modo
  cfg
}

# Ruta de salida dentro de salidas/<modo>/, creando la carpeta.
salida <- function(cfg, ...) {
  ruta <- file.path("salidas", cfg$modo, ...)
  dir.create(dirname(ruta), recursive = TRUE, showWarnings = FALSE)
  ruta
}

surtgis_bin <- function() {
  bin <- Sys.getenv("SURTGIS_BIN", Sys.which("surtgis"))
  if (!nzchar(bin)) stop("No se encontró el binario surtgis (define SURTGIS_BIN).")
  bin
}

# Ejecuta `surtgis --compress --output-dtype f32 <args>` (requiere SurtGIS >= 1.5) y falla
# con el mensaje de SurtGIS si hay error. Corre bajo `prlimit --as` (SURTGIS_MAX_MEM, default
# 6 GB) como resguardo: si algo se excede, falla en vez de congelar el equipo.
surtgis <- function(...) {
  args <- c("--compress", "--output-dtype", "f32", ...)
  limite <- Sys.getenv("SURTGIS_MAX_MEM", "6G")
  bytes <- as.numeric(sub("[GgMm]$", "", limite)) *
    if (grepl("[Gg]$", limite)) 1024^3 else if (grepl("[Mm]$", limite)) 1024^2 else 1
  res <- suppressWarnings(system2("prlimit", c(sprintf("--as=%.0f", bytes), shQuote(surtgis_bin()),
                                               shQuote(args)), stdout = TRUE, stderr = TRUE))
  estado <- attr(res, "status")
  if (!is.null(estado) && estado != 0) {
    stop("surtgis ", paste(args, collapse = " "), "\n", paste(res, collapse = "\n"), call. = FALSE)
  }
  invisible(res)
}

# Devuelve `ruta` tras verificar que existe (útil como último paso de targets de formato "file").
existe <- function(ruta) {
  if (!all(file.exists(ruta))) stop("No se generó: ", ruta[!file.exists(ruta)])
  ruta
}

# Elige un archivo de un target de formato "file" por su nombre base (targets no conserva los
# nombres del vector).
archivo <- function(x, nombre) {
  r <- x[basename(x) == nombre]
  if (length(r) != 1) stop("No se encontró ", nombre, " en: ", paste(basename(x), collapse = ", "))
  r
}
