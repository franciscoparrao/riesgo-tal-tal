# Fig. 1 — Study area. (a) Regional setting: Taltal and the five Atacama communes used to
# calibrate the reach model (March 2015 event); (b) Taltal commune: susceptibility coverage
# (coastal watersheds) and populated census units; (c) Taltal city: 2024 census blocks by
# population, area affected in 2015 and road network.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(terra); library(sf)
  library(tidyterra); library(ggspatial) })
source("paper/R/theme_paper.R")
setup_paper_theme("elsevier")
terraOptions(progress = 0)
G <- targets::tar_read(censo_in, store = st_corr)
leer <- function(q) { x <- st_transform(st_make_valid(st_read(G, query = q, quiet = TRUE)), 32719); st_geometry(x) <- "geometry"; x }

# (a) regional
regiones <- leer("SELECT REGION, SHAPE FROM Regional_CPV24 WHERE COD_REGION IN (2, 3)")
cuts <- c(2104, 3201, 3202, 3103, 3101, 3302)
comunas <- leer(sprintf("SELECT CUT, COMUNA, SHAPE FROM Comunal_CPV24 WHERE CUT IN (%s)", paste(cuts, collapse = ",")))
comunas$rol <- ifelse(comunas$CUT == 2104, "Taltal (independent validation)", "Calibration communes (2015 event)")
etq <- st_point_on_surface(comunas)
etq$nombre <- c("TALTAL" = "Taltal", "CHAÑARAL" = "Chañaral", "DIEGO DE ALMAGRO" = "Diego de Almagro",
                "TIERRA AMARILLA" = "Tierra Amarilla", "COPIAPÓ" = "Copiapó", "COPIAPO" = "Copiapó",
                "ALTO DEL CARMEN" = "Alto del Carmen")[toupper(etq$COMUNA)]
pa <- ggplot() +
  geom_sf(data = regiones, fill = "grey95", colour = "grey55", linewidth = 0.25) +
  geom_sf(data = comunas, aes(fill = rol), colour = "grey20", linewidth = 0.2) +
  ggrepel::geom_text_repel(data = etq, aes(label = nombre, geometry = geometry), stat = "sf_coordinates",
                           size = 2.1, min.segment.length = 0, segment.size = 0.2, box.padding = 0.25,
                           bg.color = "white", bg.r = 0.15,
                           family = "Liberation Sans") +
  scale_fill_manual(values = c("Taltal (independent validation)" = "#0072B2",
                               "Calibration communes (2015 event)" = "#56B4E9"), name = NULL) +
  annotation_scale(location = "br", height = unit(3, "pt"), text_cex = 0.55, line_width = 0.3,
                   pad_y = unit(0.9, "cm")) +
  annotation_north_arrow(location = "tr", height = unit(16, "pt"), width = unit(11, "pt"),
                         style = north_arrow_orienteering(text_size = 5)) +
  coord_sf(crs = 32719, datum = 4326) +
  scale_x_continuous(breaks = c(-71, -70, -69)) + scale_y_continuous(breaks = c(-29, -27, -25)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5), axis.line = element_blank(),
        legend.key.size = unit(7, "pt"), axis.text = element_text(size = 5.5), axis.ticks = element_line(linewidth = 0.2),
        panel.border = element_rect(fill = NA, colour = "grey50", linewidth = 0.2))
# recuadro: Chile continental con el área de (a)
chile <- st_read(G, query = "SELECT SHAPE FROM Regional_CPV24", quiet = TRUE) |> st_make_valid() |>
  st_transform(4326) |> st_crop(c(xmin = -76, xmax = -66, ymin = -56, ymax = -17)) |> st_union() |>
  st_simplify(dTolerance = 2000)
caja_a <- st_transform(st_as_sfc(st_bbox(regiones)), 4326)
inset <- ggplot() + geom_sf(data = chile, fill = "grey85", colour = NA) +
  geom_sf(data = caja_a, fill = NA, colour = "#D55E00", linewidth = 0.4) +
  coord_sf(datum = NA, expand = FALSE) + theme_void() +
  theme(panel.border = element_rect(fill = "white", colour = "grey50", linewidth = 0.2),
        panel.background = element_rect(fill = "white", colour = NA))
pa <- pa + inset_element(inset, left = 0.0, bottom = 0.52, right = 0.24, top = 1.0, align_to = "panel", ignore_tag = TRUE)

# (b) comuna de Taltal
com <- comunas[comunas$CUT == 2104, ]
susc <- rast(targets::tar_read(susc_ensamble, store = st_corr))
cobertura <- st_transform(st_as_sf(as.polygons(!is.na(aggregate(susc, 10)), dissolve = TRUE)), 32719)
cobertura <- st_intersection(cobertura[cobertura[[1]] == 1, ], st_geometry(com))
u <- st_read(targets::tar_read(unidades, store = st_corr), quiet = TRUE)
u <- u[u$poblacion > 0, ]
caja_c <- st_as_sfc(st_bbox(c(xmin = 349300, ymin = 7184900, xmax = 354100, ymax = 7190400), crs = 32719))
# unidades pobladas como puntos (las entidades rurales son polígonos enormes y engañan)
up <- st_point_on_surface(u[, "poblacion"]); st_geometry(up) <- "geometry"
up <- rbind(up[u$poblacion > 0 & as.numeric(st_area(u)) > 5e4, ],
            st_sf(poblacion = sum(u$poblacion[as.numeric(st_area(u)) <= 5e4]),
                  geometry = st_centroid(st_union(u[as.numeric(st_area(u)) <= 5e4, ]))))
pb <- ggplot() +
  geom_sf(data = com, fill = "grey97", colour = "grey30", linewidth = 0.3) +
  geom_sf(data = cobertura, aes(fill = "Extent of the initiation\nsusceptibility map"), alpha = 0.3, colour = NA) +
  geom_sf(data = up, aes(size = poblacion), shape = 21, fill = "#E69F00", colour = "#8A5500", stroke = 0.3) +
  geom_sf(data = caja_c, fill = NA, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = "#CC79A7", name = NULL) +
  scale_size_area(max_size = 5, breaks = c(10, 100, 1000, 10000), name = "Residents (2024)") +
  annotation_scale(location = "bl", height = unit(3, "pt"), text_cex = 0.55, line_width = 0.3) +
  annotation_north_arrow(location = "tr", height = unit(16, "pt"), width = unit(11, "pt"),
                         style = north_arrow_orienteering(text_size = 5)) +
  coord_sf(crs = 32719, datum = NA) +
  guides(fill = guide_legend(order = 1), size = guide_legend(order = 2, nrow = 1)) +
  theme(axis.line = element_blank(), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        legend.spacing.y = unit(1, "pt"), legend.key.size = unit(7, "pt"))

# (c) ciudad
caja <- ext(349300, 354100, 7184900, 7190400)
dem <- crop(rast("salidas/corregido/amenaza/dem.tif"), caja)
hs <- shade(terrain(dem, "slope", unit = "radians"), terrain(dem, "aspect", unit = "radians"), 40, 315)
names(hs) <- "hs"
mz <- leer("SELECT n_per, SHAPE FROM Manzanas_CPV24 WHERE CUT = 2104 AND n_per > 0")
afect <- st_read("salidas/corregido/validacion/area_afectada_2015.gpkg", quiet = TRUE)
vias <- st_transform(st_read("salidas/corregido/red/taltal.osm.pbf", quiet = TRUE,
                             query = "SELECT highway FROM lines WHERE highway IS NOT NULL"), 32719)
vias <- st_crop(vias, st_bbox(caja_c))
pc <- ggplot() +
  geom_spatraster(data = hs, show.legend = FALSE) +
  scale_fill_gradient(low = "grey40", high = "grey97", na.value = NA, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = mz, aes(fill = n_per), colour = NA, alpha = 0.9) +
  geom_sf(data = vias, colour = "grey25", linewidth = 0.12) +
  geom_sf(data = afect, aes(colour = "Area affected, March 2015"), fill = "#0072B2", alpha = 0.25, linewidth = 0.35) +
  scale_colour_manual(values = "#0072B2", name = NULL) +
  scico::scale_fill_scico(palette = "lajolla", direction = -1, trans = "sqrt", breaks = c(25, 100, 300, 600),
                          name = "Residents per block (2024)",
                          guide = guide_colourbar(barwidth = unit(55, "pt"), barheight = unit(5, "pt"),
                                                  title.position = "top")) +
  annotation_scale(location = "bl", height = unit(3, "pt"), text_cex = 0.55, line_width = 0.3) +
  annotation_north_arrow(location = "tr", height = unit(16, "pt"), width = unit(11, "pt"),
                         style = north_arrow_orienteering(text_size = 5)) +
  coord_sf(crs = 32719, xlim = c(xmin(caja), xmax(caja)), ylim = c(ymin(caja), ymax(caja)),
           expand = FALSE, datum = NA) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(55, "pt"), barheight = unit(5, "pt"),
                                title.position = "top"), colour = guide_legend(order = 2)) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 7),
        legend.text = element_text(size = 6), legend.spacing.y = unit(1, "pt"),
        axis.line = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3))

etq_panel <- function(p, l) p + labs(tag = l) + theme(plot.tag = element_text(face = "bold", size = 10))
fig <- (etq_panel(pa, "a)") | etq_panel(pb, "b)") | etq_panel(pc, "c)")) + plot_layout(widths = c(1, 1, 1.05))
save_fig(fig, "fig1_area_estudio", width_cm = 18, height_cm = 10)
