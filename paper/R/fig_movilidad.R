# Fig. 2 (versión movilidad) — The March 2015 aluviones and the Flow-R friction model.
# (a) Channel gradient (500 m downstream, D8, Copernicus GLO-30) of the channel cells inside the
#     affected areas of the six communes (descriptive).
# (b) Critical travel angle of the affected census blocks under the Flow-R simplified friction-limited
#     model propagated from Flow-R source areas over the whole catchments (largest angle at which a
#     flow reaches the block), with the 15 m/s velocity limit and without it (pure energy line).
#     Labels: % of affected blocks not reached at any angle >= 0.25°.
# Also writes the supplementary histogram (Fig. S2) of the travel angle needed by the susceptible
# sources of the Quebrada de Taltal to reach the mapped area.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("paper/R/theme_paper.R")
setup_paper_theme("elsevier")

nombres <- c("TALTAL" = "Taltal", "CHAÑARAL" = "Chañaral", "DIEGO DE ALMAGRO" = "Diego de Almagro",
             "TIERRA AMARILLA" = "Tierra Amarilla", "COPIAPO" = "Copiapó", "ALTO DEL CARMEN" = "Alto del Carmen")
orden <- c("TALTAL", "CHAÑARAL", "DIEGO DE ALMAGRO", "TIERRA AMARILLA", "COPIAPO", "ALTO DEL CARMEN")
niv <- rev(nombres[orden])
umbrales <- data.frame(x = c(7, 11), lab = c("7° fine", "11° coarse"))
lineas <- list(
  geom_vline(data = umbrales, aes(xintercept = x), linetype = "22", linewidth = 0.35, colour = "grey30"),
  scale_x_continuous(limits = c(0, 17.6), breaks = seq(0, 15, 3), expand = expansion(mult = c(0, 0.01)))
)

g <- read.csv("salidas/corregido/movilidad/pendiente_celdas.csv")
r <- read.csv("salidas/corregido/movilidad/pendiente_cauces_afectados.csv")
km <- setNames(ifelse(r$largo_cauce_km < 10, sprintf("%.1f", r$largo_cauce_km), sprintf("%.0f", r$largo_cauce_km)), r$comuna)
g$comuna <- factor(nombres[g$comuna], niv)
pa <- ggplot(g, aes(gradiente, comuna)) +
  geom_violin(fill = "#56B4E9", colour = NA, alpha = 0.55, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  annotate("text", x = 17.5, y = seq_along(niv), label = paste(km[names(niv)], "km"), hjust = 1, size = 2.3,
           colour = "grey30") +
  lineas +
  annotate("text", x = umbrales$x - 0.15, y = 6.45, label = paste("Flow-R", umbrales$lab), hjust = 1,
           size = 2.2, colour = "grey30") +
  labs(x = "Channel gradient in the affected reaches (°)", y = NULL) +
  coord_cartesian(clip = "off")

d <- readRDS("salidas/corregido/movilidad/bloques_sflm.rds")
d <- d[d$afectada, ]
l <- rbind(data.frame(comuna = d$comuna, modelo = "With velocity limit (15 m/s)", angulo = d$a15),
           data.frame(comuna = d$comuna, modelo = "Without limit (energy line, H/L)", angulo = d$ainf))
l$comuna <- factor(nombres[l$comuna], niv)
no <- aggregate(angulo ~ comuna + modelo, l, function(v) round(100 * mean(v == 0)))
no$lab <- paste0(no$angulo, "%")
pb <- ggplot(l[l$angulo > 0, ], aes(angulo, comuna, colour = modelo)) +
  geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.35, position = position_dodge(0.7)) +
  geom_text(data = no, aes(x = 17.5, y = comuna, label = lab, group = modelo), hjust = 1, size = 2.3,
            position = position_dodge(0.7), show.legend = FALSE) +
  lineas +
  annotate("text", x = 17.5, y = 6.62, label = "not\nreached", hjust = 1, vjust = 0, size = 2.2, lineheight = 0.85,
           colour = "grey30") +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = c("#0072B2", "#E69F00"), name = NULL) +
  labs(x = "Critical travel angle of the affected blocks (°)", y = NULL) +
  theme(legend.position = "top")

fig <- pa / pb + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))
save_fig(fig, "fig_movilidad", width_cm = 14, height_cm = 16)

a <- read.csv("salidas/corregido/flowpy/angulo_requerido_fuentes.csv")
ps <- ggplot(a, aes(angulo)) +
  geom_histogram(binwidth = 0.25, boundary = 0, fill = "#E69F00", colour = "white", linewidth = 0.1) +
  lineas +
  labs(x = "Travel angle needed to reach the mapped area of Taltal from the susceptible sources (°)",
       y = "Source cells")
save_fig(ps, "figS2_angulo_requerido_taltal", width_cm = 14, height_cm = 6.5)
