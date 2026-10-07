# Regenera todas las figuras del paper desde la raíz del proyecto (requiere el pipeline
# `Rscript run.R corregido` ya corrido):  Rscript paper/R/todas_las_figuras.R
for (f in sprintf("paper/R/fig%d_%s.R", 1:5,
                  c("area_estudio", "validacion_taltal", "validacion_externa", "mapas_riesgo", "sensibilidad"))) {
  message("== ", f)
  local(source(f, local = TRUE))
}
sessionInfo()
