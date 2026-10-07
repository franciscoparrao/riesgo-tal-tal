# Fig. 4 — Taltal city and ravine: (a) reach index; (b) reach-based multiplicative risk
# (this study); (c) conventional additive index. Outline: area affected in March 2015.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(terra); library(sf)
  library(tidyterra); library(ggspatial) })
source("paper/R/theme_paper.R")
setup_paper_theme("elsevier")
terraOptions(progress = 0)

marco <- st_read("salidas/corregido/validacion/marco_2015.gpkg", quiet = TRUE)
caja <- ext(349300, 354100, 7184900, 7190400)
afect <- st_read("salidas/corregido/validacion/area_afectada_2015.gpkg", quiet = TRUE)
# zona de riesgo por aluvión (ARA) del Plan Regulador Comunal de Taltal (2012)
ara <- st_transform(st_make_valid(st_read("datos/entrada/minvu/prc_riesgo_taltal.geojson", quiet = TRUE)), 32719)
ara <- st_union(ara[ara$ZONA == "ARA", ])
C <- "salidas/corregido/componentes/"; R <- "salidas/corregido/riesgo/"

leer <- function(f) crop(rast(f), caja)
dem <- crop(rast("salidas/corregido/amenaza/dem.tif"), caja)
hs <- shade(terrain(dem, "slope", unit = "radians"), terrain(dem, "aspect", unit = "radians"), 40, 315)
names(hs) <- "hs"
palc <- leer(paste0(C, "p_alcance.tif"))
# el riesgo solo está definido en unidades pobladas; P(alcance) en todo el terreno
names(palc) <- "v"
# solo continente: en el mar HAND = 0 por la corrección costera (P ~ 0.9 sin sentido físico)
tierra <- resample(dem, palc) > 0
palc <- mask(palc, tierra, maskvalues = c(0, NA))
# (b, c): valores por unidad censal (la unidad de análisis), como percentil entre las 154 unidades
# de validación, para que ambos índices se comparen en la misma escala; solo manzanas urbanas.
d <- targets::tar_read(tabla_2015, store = st_corr)
C3 <- (d$v_extendida_indice + d$exposicion + d$locc) / 3
d$riesgo <- d$p_alcance * C3
d$aditivo <- 0.25 * (d$susceptibilidad + d$v_extendida_indice + d$exposicion + d$locc)
pct <- function(x) 100 * (rank(x) - 0.5) / length(x)
d$p_riesgo <- pct(d$riesgo); d$p_aditivo <- pct(d$aditivo)
u <- st_read(targets::tar_read(unidades, store = st_corr), quiet = TRUE)
u <- merge(u[as.numeric(st_area(u)) < 1e6, "unidad"], d[, c("unidad", "p_riesgo", "p_aditivo")], by = "unidad")
base <- function(r, pal, titulo, lim = NULL, dir = 1) {
  ggplot() +
    geom_spatraster(data = hs, show.legend = FALSE) +
    scale_fill_gradient(low = "grey35", high = "grey97", na.value = NA, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_spatraster(data = r, maxcell = 3e6, alpha = 0.85) +
    scico::scale_fill_scico(palette = pal, direction = dir, limits = lim, na.value = NA,
                            name = titulo, guide = guide_colourbar(barwidth = unit(55, "pt"),
                                                                   barheight = unit(5, "pt"),
                                                                   title.position = "top")) +
    geom_sf(data = afect, fill = NA, colour = "black", linewidth = 0.35, linetype = "22") +
    coord_sf(crs = 32719, xlim = c(xmin(caja), xmax(caja)), ylim = c(ymin(caja), ymax(caja)),
             expand = FALSE, datum = NA) +
    theme(legend.position = "bottom", legend.title = element_text(size = 7),
          legend.text = element_text(size = 6), axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3))
}

pa <- base(palc, "lajolla", "Reach index", c(0, 1), dir = -1) +
  geom_sf(data = st_sf(geometry = ara), aes(linetype = "Official aluvión zone ARA (2012)"), fill = NA,
          colour = "#0072B2", linewidth = 0.45) +
  scale_linetype_manual(values = "solid", name = NULL) +
  coord_sf(crs = 32719, xlim = c(xmin(caja), xmax(caja)), ylim = c(ymin(caja), ymax(caja)),
           expand = FALSE, datum = NA) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(55, "pt"), barheight = unit(5, "pt"),
                                title.position = "top"),
         linetype = guide_legend(order = 2, override.aes = list(colour = "#0072B2"))) +
  theme(legend.box = "vertical", legend.spacing.y = unit(1, "pt")) +
  annotation_scale(location = "bl", height = unit(3, "pt"), text_cex = 0.6, line_width = 0.3,
                   pad_y = unit(6, "pt")) +
  annotation_north_arrow(location = "tr", height = unit(22, "pt"), width = unit(14, "pt"),
                         style = north_arrow_orienteering(text_size = 6, fill = c("black", "white")))
base_u <- function(var, titulo) {
  ggplot() +
    geom_spatraster(data = hs, show.legend = FALSE) +
    scale_fill_gradient(low = "grey35", high = "grey97", na.value = NA, guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_sf(data = u, aes(fill = .data[[var]]), colour = "white", linewidth = 0.05, alpha = 0.9) +
    scico::scale_fill_scico(palette = "lajolla", direction = -1, limits = c(0, 100), breaks = c(0, 50, 100),
                            name = titulo, guide = guide_colourbar(barwidth = unit(55, "pt"),
                                                                   barheight = unit(5, "pt"),
                                                                   title.position = "top")) +
    geom_sf(data = afect, fill = NA, colour = "black", linewidth = 0.35, linetype = "22") +
    annotation_scale(location = "bl", height = unit(3, "pt"), text_cex = 0.6, line_width = 0.3,
                     pad_y = unit(6, "pt")) +
    annotation_north_arrow(location = "tr", height = unit(22, "pt"), width = unit(14, "pt"),
                           style = north_arrow_orienteering(text_size = 6, fill = c("black", "white"))) +
    coord_sf(crs = 32719, xlim = c(xmin(caja), xmax(caja)), ylim = c(ymin(caja), ymax(caja)),
             expand = FALSE, datum = NA) +
    theme(legend.position = "bottom", legend.title = element_text(size = 7),
          legend.text = element_text(size = 6), axis.line = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3))
}
pb <- base_u("p_riesgo", "Reach-based risk (percentile)")
pc <- base_u("p_aditivo", "Conventional additive index (percentile)")

fig <- (pa | pb | pc) + plot_annotation(tag_levels = "a", tag_suffix = ")")
save_fig(fig, "fig4_mapas_riesgo", width_cm = 18, height_cm = 9.5)
