# Fig. 2 — Initiation vs. reach in the Quebrada de Taltal (real data).
# (a) Longitudinal profile of the main channel down to the sea: susceptibility of the valley sides
#     vs. the channel floor, and the reach affected in March 2015.
# (b) Cross-section of the city perpendicular to the ravine: elevation, initiation susceptibility
#     and HAND, with the affected interval shaded.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(terra); library(sf) })
source("paper/R/theme_paper.R")
setup_paper_theme("elsevier")
terraOptions(progress = 0)

A <- "salidas/corregido/amenaza/"
entrada <- c(353917.3, 7185196)          # entrada del flujo al mapa (R/angulo_requerido.R)
caja <- ext(340000, entrada[1] + 40000, entrada[2] - 25000, entrada[2] + 25000)
rd <- function(f) crop(rast(paste0(A, f)), caja)
fd <- rd("flowdir_terra.tif"); dem <- rd("dem_fill.tif"); s <- rd("susc_30m.tif")
hand <- rd("hand_costa_2000.tif")
afect <- st_union(st_read("salidas/corregido/validacion/area_afectada_2015.gpkg", quiet = TRUE))

# --- red D8 (terra flowdir: 1=E 2=SE 4=S 8=SW 16=W 32=NW 64=N 128=NE) ---
nc <- ncol(fd); nr <- nrow(fd)
k <- match(values(fd, mat = FALSE), c(1, 2, 4, 8, 16, 32, 64, 128))
dc <- c(1, 1, 0, -1, -1, -1, 0, 1); dr <- c(0, 1, 1, 1, 0, -1, -1, -1)
paso <- 30 * c(1, sqrt(2), 1, sqrt(2), 1, sqrt(2), 1, sqrt(2))
i0 <- seq_along(k) - 1L; r0 <- i0 %/% nc; c0 <- i0 %% nc
r1 <- r0 + dr[k]; c1 <- c0 + dc[k]
abajo <- ifelse(!is.na(k) & r1 >= 0 & r1 < nr & c1 >= 0 & c1 < nc, r1 * nc + c1 + 1L, NA_integer_)
acc <- values(flowAccumulation(fd), mat = FALSE)
z <- values(dem, mat = FALSE)

# cauce principal: aguas arriba desde la entrada siguiendo la mayor acumulación (35 km) y aguas
# abajo hasta el mar
ce <- cellFromXY(fd, rbind(entrada))
arriba <- ce; l <- 0
repeat {
  j <- arriba[1]; rc <- rowColFromCell(fd, j)
  vec <- cbind(rc[1] + c(-1, -1, -1, 0, 0, 1, 1, 1), rc[2] + c(-1, 0, 1, -1, 1, -1, 0, 1))
  vec <- vec[vec[, 1] >= 1 & vec[, 1] <= nr & vec[, 2] >= 1 & vec[, 2] <= nc, , drop = FALSE]
  v <- cellFromRowCol(fd, vec[, 1], vec[, 2]); v <- v[which(abajo[v] == j)]
  if (!length(v) || l > 35000) break
  u <- v[which.max(acc[v])]; l <- l + paso[k[u]]; arriba <- c(u, arriba)
}
abajo_ruta <- c(); j <- abajo[ce]
while (!is.na(j) && !is.na(z[j]) && z[j] > 1) { abajo_ruta <- c(abajo_ruta, j); j <- abajo[j] }
ruta <- c(arriba, abajo_ruta)
d_ruta <- c(0, cumsum(paso[k[ruta[-length(ruta)]]]))
costa <- max(d_ruta)
xy <- xyFromCell(fd, ruta)

# susceptibilidad de las laderas: media en un círculo de 500 m excluyendo el fondo del valle
# (HAND < 10 m); susceptibilidad del fondo: celdas del cauce
s_lad <- mask(s, hand < 10, maskvalues = 1)
w <- focalMat(s, 500, "circle"); w[w > 0] <- 1
s_lad_f <- focal(s_lad, w, fun = "mean", na.rm = TRUE)
s_cor_f <- focal(s, w, fun = "mean", na.rm = TRUE)
pa <- data.frame(dist_costa = (costa - d_ruta) / 1000, z = z[ruta],
                 laderas = extract(s_lad_f, xy)[, 1], corredor = extract(s_cor_f, xy)[, 1], cauce = values(s, mat = FALSE)[ruta],
                 afectada = lengths(st_intersects(st_as_sf(as.data.frame(xy), coords = 1:2, crs = 32719),
                                                  afect)) > 0)
pa$laderas_suav <- stats::filter(pa$laderas, rep(1 / 15, 15), sides = 2)
pa$cauce_suav <- stats::filter(pa$cauce, rep(1 / 15, 15), sides = 2)
pa$corredor_suav <- stats::filter(pa$corredor, rep(1 / 15, 15), sides = 2)
tramo <- range(pa$dist_costa[pa$afectada])

# --- (b) corte transversal en la ciudad, perpendicular al cauce, ~1.5 km de la costa ---
k0 <- which.min(abs(pa$dist_costa - 1.5))
dirv <- xy[min(k0 + 10, nrow(xy)), ] - xy[max(k0 - 10, 1), ]; dirv <- dirv / sqrt(sum(dirv^2))
perp <- c(-dirv[2], dirv[1])
t_m <- seq(-1500, 1500, by = 15)
pts <- cbind(xy[k0, 1] + t_m * perp[1], xy[k0, 2] + t_m * perp[2])
dem_c <- crop(rast(paste0(A, "dem.tif")), caja)
pb <- data.frame(t = t_m / 1000, z = extract(dem_c, pts)[, 1], s = extract(s, pts)[, 1],
                 hand = extract(hand, pts)[, 1],
                 afectada = lengths(st_intersects(st_as_sf(as.data.frame(pts), coords = 1:2, crs = 32719),
                                                  afect)) > 0)
if (pb$z[1] < pb$z[nrow(pb)]) pb$t <- -pb$t   # orientar: ladera alta a la izquierda no importa; mantener norte-sur
seg <- with(pb, { r <- rle(afectada); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  data.frame(x0 = t[ini][r$values], x1 = t[fin][r$values]) })

col_lad <- "#D55E00"; col_cau <- "#0072B2"; col_af <- "grey85"

# (a) dos filas: perfil de elevación y susceptibilidad
a1 <- ggplot(pa, aes(dist_costa, z)) +
  annotate("rect", xmin = tramo[1], xmax = tramo[2], ymin = -Inf, ymax = Inf, fill = col_af) +
  geom_line(linewidth = 0.5) +
  annotate("text", x = mean(tramo), y = Inf, vjust = 1.4, label = "affected\n2015", size = 2.5, lineheight = 0.9) +
  scale_x_reverse(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = NULL, y = "Channel elevation (m)") +
  theme(axis.text.x = element_blank())
a2 <- ggplot(pa, aes(dist_costa, corredor_suav)) +
  annotate("rect", xmin = tramo[1], xmax = tramo[2], ymin = -Inf, ymax = Inf, fill = col_af) +
  geom_line(colour = col_lad, linewidth = 0.5, na.rm = TRUE) +
  scale_x_reverse(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Distance along the channel to the coast (km)", y = "Initiation\nsusceptibility")

# (b) tres filas compartiendo x
sombra <- function() geom_rect(data = seg, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf),
                               inherit.aes = FALSE, fill = col_af)
lado <- if (perp[2] > 0) c("S", "N") else c("N", "S")   # t < 0 y t > 0
b1 <- ggplot(pb, aes(t, z)) + sombra() + geom_line(linewidth = 0.5, na.rm = TRUE) +
  annotate("text", x = mean(unlist(seg[1, ])), y = Inf, vjust = 1.4, label = "affected\n2015", size = 2.5, lineheight = 0.9) +
  labs(x = NULL, y = "Elevation (m)") + theme(axis.text.x = element_blank())
b2 <- ggplot(pb, aes(t, s)) + sombra() + geom_line(colour = col_lad, linewidth = 0.5, na.rm = TRUE) +
  annotate("text", x = c(-1.45, 1.45), y = 0.08, label = lado, size = 2.6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Initiation\nsusceptibility") + theme(axis.text.x = element_blank())
b3 <- ggplot(pb, aes(t, hand)) + sombra() + geom_line(colour = col_cau, linewidth = 0.5, na.rm = TRUE) +
  labs(x = "Distance across the ravine (km)", y = "HAND (m)")

p <- (a1 / a2 + plot_layout(heights = c(1.4, 1))) | (b1 / b2 / b3)
p <- p + plot_annotation(tag_levels = list(c("a", "", "b", "", ""))) &
  theme(plot.tag = element_text(face = "bold"))
save_fig(p, "fig_concepto", width_cm = 18, height_cm = 10)

cat(sprintf("Perfil: %.1f km; tramo afectado %.1f-%.1f km de la costa; corte a %.2f km de la costa\n",
            costa / 1000, tramo[1], tramo[2], pa$dist_costa[k0]))
cat(sprintf("Corredor 500 m: media fuera del tramo afectado %.2f, en el tramo %.2f; máximo %.2f a %.1f km de la costa\n",
            mean(pa$corredor[!pa$afectada], na.rm = TRUE), mean(pa$corredor[pa$afectada], na.rm = TRUE),
            max(pa$corredor_suav, na.rm = TRUE), pa$dist_costa[which.max(pa$corredor_suav)]))
cat(sprintf("Susceptibilidad media laderas %.2f vs cauce %.2f (todo el perfil); en el tramo afectado: %.2f vs %.2f\n",
            mean(pa$laderas, na.rm = TRUE), mean(pa$cauce, na.rm = TRUE),
            mean(pa$laderas[pa$afectada], na.rm = TRUE), mean(pa$cauce[pa$afectada], na.rm = TRUE)))
cat(sprintf("Corte: susceptibilidad media afectado %.2f vs no afectado %.2f; HAND mediana %.1f vs %.1f m\n",
            mean(pb$s[pb$afectada], na.rm = TRUE), mean(pb$s[!pb$afectada], na.rm = TRUE),
            median(pb$hand[pb$afectada], na.rm = TRUE), median(pb$hand[!pb$afectada], na.rm = TRUE)))
