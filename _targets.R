# Pipeline de riesgo por aluviones — Taltal: índice convencional, término de alcance y
# validación contra el evento de marzo de 2015.
# Ejecutar con:  Rscript run.R corregido
#
# Dependencias de configuración: ningún target costoso menciona `cfg` en su comando (targets
# lo trataría como dependencia del objeto completo y cualquier cambio en config.yml lo
# invalidaría). Reciben sub-configuraciones `p_*`, que son targets baratos: si un campo cambia,
# solo se invalida lo que usa ese `p_*`. Todas incluyen `modo` y `crs` (salida() usa `modo`).

library(targets)
library(tarchetypes)

tar_source("R")
tar_option_set(packages = c("terra", "sf", "classInt", "yaml", "data.table", "dbscan", "cluster", "aricode", "e1071"),
               format = "rds")

cfg <- leer_config()
parte <- function(...) c(cfg[c("modo", "crs")], cfg[c(...)])

comunes <- list(
  tar_target(config, cfg),
  tar_target(p_base, parte()),
  tar_target(p_grilla, parte("resolucion_salida", "resolucion_locc")),
  tar_target(p_integracion, parte("pesos", "semilla", "jenks_clases", "jenks_muestra")),

  # --- Entradas (rutas: se re-evalúan baratas; si el archivo no cambia, nada aguas abajo) ---
  tar_target(comuna_in, cfg$entradas$comuna, format = "file")
)

# --- Modo corregido ------------------------------------------------------------------------
indicadores <- data.frame(
  indicador = names(cfg$locc_indicadores),
  capa = vapply(cfg$locc_indicadores, `[[`, "", "capa"),
  umbral = vapply(cfg$locc_indicadores, `[[`, 0, "umbral")
)

combinaciones <- expand.grid(variante = cfg$vulnerabilidad$variantes,
                             metodo = cfg$vulnerabilidad$metodos, stringsAsFactors = FALSE)
tipologias <- expand.grid(variante = cfg$vulnerabilidad$variantes,
                          metodo = cfg$vulnerabilidad$tipologias, stringsAsFactors = FALSE)

modo_corregido <- list(
  tar_target(p_red, parte("osm_margen_km")),
  tar_target(p_vulnerabilidad, parte("vulnerabilidad", "semilla", "censo")),
  tar_target(p_exposicion, parte("censo", "exposicion")),

  # Red vial y grilla fina
  tar_target(osm_in, cfg$osm_pbf, format = "file"),
  tar_target(osm_taltal, recortar_osm(osm_in, comuna_in, p_red$osm_margen_km, p_red), format = "file"),
  tar_target(red, red_r5r(osm_taltal), format = "file"),
  tar_target(recreativos, recreativos_osm(osm_taltal, comuna_in, p_base), format = "file"),
  tar_target(comuna_utm, {
    d <- salida(p_base, "tmp", "comuna_32719.gpkg")
    sf::st_write(sf::st_transform(sf::st_read(comuna_in, quiet = TRUE), p_base$crs), d,
                 delete_dsn = TRUE, quiet = TRUE)
    d
  }, format = "file"),
  tar_target(plantilla_fina,
             plantilla(comuna_utm, p_grilla$resolucion_locc, salida(p_base, "tmp", "plantilla_10m.tif")),
             format = "file"),

  # LOCC
  tar_map(
    values = indicadores, names = "indicador",
    tar_target(capa_in, if (capa == "osm") recreativos else capa, format = "file"),
    tar_target(instalaciones, leer_instalaciones(capa_in)),
    tar_target(isocrona, isocrona_tipo(red, instalaciones, umbral, indicador, p_base), format = "file"),
    tar_target(indicador_locc_r, indicador_locc(isocrona, comuna_in, plantilla_fina, indicador, p_base),
               format = "file")
  ),
  tar_target(locc, locc_promedio(c(indicador_locc_r_salud, indicador_locc_r_comisarias,
                                   indicador_locc_r_educacion, indicador_locc_r_transmision,
                                   indicador_locc_r_recreacion), p_base), format = "file"),

  # Exposición
  tar_target(censo_in, cfg$censo$ruta, format = "file"),
  tar_target(buildings_in, cfg$google_buildings, format = "file"),
  tar_target(exposicion_piezas, exposicion_vector(censo_in, buildings_in, osm_taltal, comuna_utm, p_exposicion),
             format = "file"),
  tar_target(exposicion, exposicion_raster(exposicion_piezas, plantilla_fina, p_base), format = "file"),

  # Vulnerabilidad
  tar_target(s2_principal, file.path(cfg$sentinel2_dir, paste0(bandas_s2, ".tif")), format = "file"),
  tar_target(s2_costa, file.path(cfg$sentinel2_dir, "relleno_costa", paste0(bandas_s2, ".tif")), format = "file"),
  tar_target(s2_bandas, rellenar_s2(s2_principal, s2_costa, p_base), format = "file"),
  tar_target(unidades, unidades_archivo(censo_in, comuna_utm, p_vulnerabilidad), format = "file"),
  tar_target(muestra, muestra_vulnerabilidad(unidades, dirname(s2_bandas[1]), p_vulnerabilidad)),
  tar_target(geologia_in, cfg$geologia, format = "file"),
  tar_target(diagnostico, diagnostico_espectral(muestra, geologia_in, p_vulnerabilidad), format = "file"),
  tar_target(espectral, espectral_unidades(unidades, dirname(s2_bandas[1]), p_vulnerabilidad)),
  # Tipologías descriptivas (K-means, HDBSCAN): solo perfiles
  tar_map(
    values = tipologias, names = c("variante", "metodo"),
    tar_target(tipologia, modelo_vulnerabilidad(unidades, espectral, variante, metodo, p_vulnerabilidad)),
    tar_target(perfiles_tipologia, perfiles_csv(tipologia, p_base), format = "file")
  ),

  # Riesgo: índice continuo vs fuzzy c-means, x variante
  tar_target(susc_modelos, unlist(cfg$susceptibilidad_modelos), format = "file"),
  tar_target(susc_ensamble, susceptibilidad_ensamble(susc_modelos, p_base), format = "file"),
  tar_target(susceptibilidad, susceptibilidad_fina(susc_ensamble, plantilla_fina, p_base), format = "file"),
  tar_target(puntos_pobladas, puntos_unidades(unidades, p_vulnerabilidad)),
  mapa_riesgo <- tar_map(
    values = combinaciones, names = c("variante", "metodo"),
    tar_target(modelo, modelo_vulnerabilidad(unidades, espectral, variante, metodo, p_vulnerabilidad)),
    tar_target(perfiles, perfiles_csv(modelo, p_base), format = "file"),
    tar_target(vulnerabilidad, vulnerabilidad_raster(modelo, unidades, plantilla_fina, p_base),
               format = "file"),
    tar_target(riesgo, riesgo_corregido(susceptibilidad, vulnerabilidad, exposicion, locc, p_integracion,
                                        paste(variante, metodo, sep = "_")), format = "file"),
    tar_target(resumen, resumen_riesgo(riesgo, puntos_pobladas, p_integracion, paste(variante, metodo, sep = "_")))
  ),
  tar_combine(comparacion, mapa_riesgo[["resumen"]],
              command = comparar_riesgos(list(!!!.x), p_base), format = "file"),

  # Validación: aluvión de marzo de 2015 (Sernageomin, digitalizado de capturas)
  tar_target(png_2015_norte, cfg$evento_2015$norte, format = "file"),
  tar_target(png_2015_sur, cfg$evento_2015$sur, format = "file"),
  tar_target(area_afectada_2015, digitalizar_2015(png_2015_norte, png_2015_sur, p_base), format = "file"),
  tar_target(marco_2015_v, marco_2015(png_2015_norte, png_2015_sur, p_base), format = "file"),
  # Validación externa (mismo evento, 6 comunas) y calibración de P(alcance) SOLO con Atacama
  tar_target(catastro_in, cfg$evento_2015$catastro, format = "file"),
  tar_target(manzanas_validacion, validacion_externa_datos(catastro_in, censo_in, area_afectada_2015, p_base)),
  tar_target(auc_comunas, guardar_csv(
    auc_tabla(manzanas_validacion, list(HAND = function(x) -x$hand, "distancia al cauce" = function(x) -x$dist),
              p_integracion$semilla),
    p_base, "validacion_externa", "auc_comunas.csv"), format = "file"),
  tar_target(comunas_calibracion, setdiff(comunas_validacion$comuna, "TALTAL")),
  tar_target(modelo_alcance, calibrar_alcance_final(manzanas_validacion, comunas_calibracion)),
  tar_target(sensibilidad_alcance_csv, guardar_csv(
    sensibilidad_alcance(manzanas_validacion, comunas_calibracion, p_integracion$semilla),
    p_base, "validacion_externa", "sensibilidad_alcance.csv"), format = "file"),

  # Amenaza por flujo (fuente -> propagación) y riesgo multiplicativo (resultado principal)
  tar_target(amenaza, amenaza_flujo(comuna_utm, susc_ensamble, modelo_alcance, plantilla_fina, p_base),
             format = "file"),
  tar_target(riesgo_final, riesgo_multiplicativo(archivo(amenaza, "p_alcance.tif"), vulnerabilidad_extendida_indice,
                                                 exposicion, locc, p_base, "final"), format = "file"),

  # Validación en Taltal 2015 (independiente de la calibración)
  tar_target(validacion_2015, validar_2015(
    area_afectada_2015, marco_2015_v, unidades, puntos_pobladas,
    capas = c(susceptibilidad = susceptibilidad, p_alcance = archivo(amenaza, "p_alcance.tif"),
              amenaza_flujo_media = archivo(amenaza, "amenaza_flujo.tif"), amenaza_fuente_total = archivo(amenaza, "amenaza_fuente_total.tif"),
              amenaza_combinada = archivo(amenaza, "amenaza_combinada.tif"),
              exposicion = exposicion, locc = locc,
              vulnerabilidad_extendida_indice = vulnerabilidad_extendida_indice,
              vulnerabilidad_base_indice = vulnerabilidad_base_indice,
              riesgo_aditivo_tesis = riesgo_extendida_indice,
              riesgo_final = riesgo_final),
    p = p_integracion), format = "file"),
  # Propagación (susceptibilidad aguas abajo, LAHARZ, cono) vs alcance, mismo protocolo
  tar_target(propagacion, capas_propagacion(amenaza, marco_2015_v, p_base), format = "file"),
  tar_target(validacion_propagacion_2015, validar_2015(
    area_afectada_2015, marco_2015_v, unidades, puntos_pobladas,
    capas = c(p_alcance = archivo(amenaza, "p_alcance.tif"), hand_neg = archivo(propagacion, "hand_neg.tif"),
              dist_cauce_neg = archivo(propagacion, "dist_cauce_neg.tif"),
              susceptibilidad = susceptibilidad,
              susc_aguas_arriba = archivo(propagacion, "susc_aguas_arriba_30m.tif"),
              fuente_susceptible = archivo(propagacion, "fuente_susceptible.tif"),
              laharz_graduado = archivo(propagacion, "laharz_graduado.tif"),
              cono_graduado = archivo(propagacion, "cono_graduado.tif")),
    p = p_integracion, nombre = "validacion_propagacion_2015"), format = "file"),
  # Flow-Py (Flow-R) sembrado desde la susceptibilidad, parámetros publicados; ángulo requerido
  tar_target(flowpy_entradas, entradas_flowpy(amenaza, marco_2015_v, p_base), format = "file"),
  tar_target(flowpy_a11, correr_flowpy(flowpy_entradas, 11), format = "file"),
  tar_target(flowpy_a7, correr_flowpy(flowpy_entradas, 7), format = "file"),
  tar_target(validacion_flowpy_2015, validar_2015(
    area_afectada_2015, marco_2015_v, unidades, puntos_pobladas,
    capas = c(p_alcance = archivo(amenaza, "p_alcance.tif"), flowpy_a11 = flowpy_a11, flowpy_a7 = flowpy_a7),
    p = p_integracion, nombre = "validacion_flowpy_2015"), format = "file"),
  tar_target(angulo_requerido_res, angulo_requerido(amenaza, marco_2015_v, p_base)),
  tar_target(angulo_requerido_csv, c(
    guardar_csv(angulo_requerido_res$resumen, p_base, "flowpy", "angulo_requerido.csv"),
    guardar_csv(angulo_requerido_res$angulos, p_base, "flowpy", "angulo_requerido_fuentes.csv")), format = "file"),
  # Movilidad: pendiente de los cauces afectados en las seis comunas
  tar_target(movilidad_res, movilidad_comunas(catastro_in, area_afectada_2015, manzanas_validacion, p_base)),
  tar_target(movilidad_csv, c(
    guardar_csv(movilidad_res$resumen, p_base, "movilidad", "pendiente_cauces_afectados.csv"),
    guardar_csv(movilidad_res$celdas, p_base, "movilidad", "pendiente_celdas.csv")), format = "file"),
  # SFLM de Flow-R en cuencas completas (bin/run_sflm.sh, ~1 h, fuera de targets) y evaluación
  tar_target(sflm_rasters, Sys.glob("salidas/corregido/movilidad/region/*/angulo_v*.tif"), format = "file"),
  tar_target(sflm_bloques, { stopifnot(length(sflm_rasters) == 8)
    bloques_sflm(catastro_in, censo_in, area_afectada_2015, "salidas/corregido/movilidad/region", p_base) }),
  tar_target(sflm_eval, evaluar_sflm(sflm_bloques, p_integracion$semilla)),
  tar_target(sflm_csv, c(guardar_csv(sflm_eval$por_comuna, p_base, "movilidad", "sflm_por_comuna.csv"),
                         guardar_csv(sflm_eval$loco, p_base, "movilidad", "sflm_loco.csv")), format = "file"),
  # Líneas base topográficas triviales (elevación, pendiente, distancia a la costa)
  tar_target(triviales, capas_triviales_taltal(amenaza, p_base), format = "file"),
  tar_target(validacion_trivial_2015, validar_2015(
    area_afectada_2015, marco_2015_v, unidades, puntos_pobladas,
    capas = c(p_alcance = archivo(amenaza, "p_alcance.tif"), hand_neg = archivo(propagacion, "hand_neg.tif"),
              elev_neg = archivo(triviales, "elev_neg.tif"), pendiente_neg = archivo(triviales, "pendiente_neg.tif"),
              dist_costa_neg = archivo(triviales, "dist_costa_neg.tif")),
    p = p_integracion, nombre = "validacion_trivial_2015"), format = "file"),
  tar_target(triviales_comunas_csv, guardar_csv(
    auc_triviales_comunas(catastro_in, censo_in, area_afectada_2015, manzanas_validacion, p_base),
    p_base, "validacion_externa", "auc_triviales.csv"), format = "file"),
  tar_target(validacion_susc_2015,
             validar_susceptibilidad_2015(area_afectada_2015, marco_2015_v, susceptibilidad, p_base),
             format = "file"),

  # Análisis de sensibilidad
  tar_target(sens_umbral_csv, guardar_csv(
    sens_umbral(catastro_in, censo_in, area_afectada_2015, comunas_calibracion, p_integracion),
    p_base, "sensibilidad", "umbral_cauce.csv"), format = "file"),
  tar_target(fabdem_in, cfg$dem_fabdem, format = "file"),
  tar_target(sens_dem_csv, guardar_csv(
    sens_dem(censo_in, area_afectada_2015, fabdem_in, modelo_alcance, manzanas_validacion, p_integracion),
    p_base, "sensibilidad", "dem.csv"), format = "file"),
  tar_target(tabla_2015, tabla_unidades_2015(
    area_afectada_2015, marco_2015_v, unidades, puntos_pobladas,
    capas = c(p_alcance = archivo(amenaza, "p_alcance.tif"), susceptibilidad = susceptibilidad,
              v_extendida_indice = vulnerabilidad_extendida_indice, v_base_indice = vulnerabilidad_base_indice,
              v_extendida_fuzzy = vulnerabilidad_extendida_fuzzy, v_base_fuzzy = vulnerabilidad_base_fuzzy,
              exposicion = exposicion, locc = locc))),
  tar_target(sens_riesgo_res, sens_riesgo(tabla_2015, p_integracion)),
  # Robustez de la validación (bootstrap espacial, correlaciones parciales, β, censo 2017)
  tar_target(robustez_res, robustez_2015(
    tabla_2015, unidades, puntos_pobladas, file.path(dirname(salida(p_base, "amenaza", "x")), "hand_costa_2000.tif"),
    manzanas_validacion, comunas_calibracion, modelo_alcance, p_integracion)),
  tar_target(robustez_csv, vapply(names(robustez_res), function(n)
    guardar_csv(robustez_res[[n]], p_base, "robustez", paste0(n, ".csv")), ""), format = "file"),
  # Robustez II: huella desplazada, nivel base marino, buffer de calibración, AUC pareados
  tar_target(sens_huella_csv, guardar_csv(sens_huella(tabla_2015, area_afectada_2015, unidades, semilla = p_integracion$semilla),
                                          p_base, "robustez", "huella_desplazada.csv"), format = "file"),
  tar_target(sens_hand_buffer_csv, guardar_csv(
    sens_hand_buffer(catastro_in, censo_in, area_afectada_2015, manzanas_validacion, p_base),
    p_base, "robustez", "hand_buffer_nivel_base.csv"), format = "file"),
  tar_target(auc_pareado_csv, guardar_csv(auc_pareado(
    tabla_pareada(c(validacion_2015, validacion_propagacion_2015, validacion_trivial_2015, validacion_flowpy_2015)),
    unidades, list(c("p_alcance", "dist_cauce_neg"), c("p_alcance", "elev_neg"), c("p_alcance", "pendiente_neg"),
                   c("p_alcance", "laharz_graduado"), c("p_alcance", "flowpy_a7"), c("p_alcance", "susceptibilidad"),
                   c("riesgo_final", "p_alcance"), c("riesgo_final", "riesgo_aditivo_tesis")),
    semilla = p_integracion$semilla), p_base, "robustez", "auc_pareado.csv"), format = "file"),
  # Robustez III (review Natural Hazards)
  tar_target(riesgo_sin_escala_csv, guardar_csv(riesgo_sin_escala(tabla_2015, p_integracion$semilla),
                                                p_base, "robustez", "riesgo_sin_escala.csv"), format = "file"),
  tar_target(bloques_coord, bloques_con_coordenadas(catastro_in, censo_in, area_afectada_2015, p_base)),
  tar_target(auc_bloques_comunas_csv, guardar_csv(auc_bloques_comunas(bloques_coord, p_integracion$semilla),
                                                  p_base, "robustez", "auc_bloques_comunas.csv"), format = "file"),
  tar_target(seleccion_bootstrap_csv, guardar_csv(seleccion_bootstrap(manzanas_validacion, comunas_calibracion,
                                                  semilla = p_integracion$semilla), p_base, "robustez", "seleccion_bootstrap.csv"), format = "file"),
  tar_target(taltal_tipo_csv, guardar_csv(taltal_por_tipo(bloques_coord, p_integracion$semilla),
                                          p_base, "robustez", "taltal_por_tipo.csv"), format = "file"),
  tar_target(holm_csv, guardar_csv(holm_pareado(auc_pareado_csv), p_base, "robustez", "auc_pareado_holm.csv"), format = "file"),
  # Selección ampliada (GFI, elevación, pendiente), contrastes con las 6 comunas, barrido aditivo
  tar_target(area_cauce_rasters, Sys.glob("salidas/corregido/movilidad/region/*/area_cauce_2000.tif"), format = "file"),
  tar_target(bloques_amp, { stopifnot(length(area_cauce_rasters) == 4); triviales_comunas_csv
    bloques_ampliado(bloques_coord, p_base) }),
  tar_target(seleccion_ampliada_csv, guardar_csv(seleccion_ampliada(bloques_amp, p_integracion$semilla),
                                                 p_base, "robustez", "seleccion_ampliada.csv"), format = "file"),
  tar_target(contrastes_res, contrastes_conjuntos(bloques_amp, semilla = p_integracion$semilla)),
  tar_target(contrastes_csv, c(guardar_csv(contrastes_res$auc_media, p_base, "robustez", "contrastes_auc_media.csv"),
                               guardar_csv(contrastes_res$diferencias, p_base, "robustez", "contrastes_diferencias.csv"),
                               guardar_csv(contrastes_res$mixtos, p_base, "robustez", "contrastes_mixtos.csv")), format = "file"),
  tar_target(barrido_aditivo_csv, guardar_csv(barrido_aditivo(tabla_2015), p_base, "robustez", "barrido_aditivo.csv"),
             format = "file"),
  # Zonificación oficial (PRC Taltal 2012, anterior al evento) como término de amenaza de planificador
  tar_target(prc_riesgo_in, "datos/entrada/minvu/prc_riesgo_taltal.geojson", format = "file"),
  tar_target(limite_urbano_in, "datos/entrada/minvu/limite_urbano_taltal.geojson", format = "file"),
  tar_target(zonificacion_res, zonificacion_oficial(prc_riesgo_in, tabla_2015, unidades, area_afectada_2015,
                                                    limite_urbano_in, semilla = p_integracion$semilla)),
  tar_target(zonificacion_csv, c(guardar_csv(zonificacion_res$validacion, p_base, "robustez", "zonificacion_oficial.csv"),
                                 guardar_csv(zonificacion_res$pareado, p_base, "robustez", "zonificacion_pareado.csv"),
                                 guardar_csv(zonificacion_res$fuera_ara, p_base, "robustez", "zonificacion_fuera_ara.csv")),
             format = "file"),
  # Paquete para Zenodo
  tar_target(deposito, exportar_deposito(tabla_2015, unidades, bloques_amp,
    file.path(dirname(salida(p_base, "amenaza", "x")), "p_alcance_30m.tif"),
    c(validacion_2015, validacion_propagacion_2015, validacion_trivial_2015, validacion_flowpy_2015,
      seleccion_ampliada_csv, contrastes_csv, zonificacion_csv, barrido_aditivo_csv, auc_pareado_csv,
      sens_umbral_csv, sens_dem_csv, sens_riesgo_csv, robustez_csv), p_base), format = "file"),
  tar_target(zonas_2017_in, "datos/entrada/censo/taltal_geo.shp", format = "file"),
  tar_target(censo_2017_2024_csv, guardar_csv(
    censo_2017_vs_2024(zonas_2017_in, unidades, area_afectada_2015, marco_2015_v, p_base),
    p_base, "robustez", "censo_2017_vs_2024.csv"), format = "file"),
  tar_target(sens_riesgo_csv, c(
    guardar_csv(sens_riesgo_res$formulaciones, p_base, "sensibilidad", "formulaciones_riesgo.csv"),
    guardar_csv(sens_riesgo_res$monte_carlo, p_base, "sensibilidad", "monte_carlo_pesos.csv")), format = "file")
)

c(comunes, modo_corregido)
