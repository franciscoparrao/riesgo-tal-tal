# Fig. 3 — External validation: same March 2015 event, six communes of the Atacama Desert.
# AUC (95 % bootstrap CI) of HAND, distance to channel and two trivial topographic baselines
# (elevation, slope) for affected vs non-affected census blocks. The reach model is calibrated on
# the five Atacama communes only; Taltal is independent. Rows: commune and geomorphic setting.
suppressPackageStartupMessages({ library(ggplot2) })
source("paper/R/theme_paper.R")
setup_paper_theme("elsevier")

t <- read.csv("salidas/corregido/validacion_externa/auc_triviales.csv", check.names = FALSE)
t <- t[t$predictor != "distance to coast", ]
nombres <- c("TALTAL" = "Taltal", "CHAÑARAL" = "Chañaral", "DIEGO DE ALMAGRO" = "Diego de Almagro",
             "TIERRA AMARILLA" = "Tierra Amarilla", "COPIAPO" = "Copiapó", "ALTO DEL CARMEN" = "Alto del Carmen")
entorno <- c("TALTAL" = "coastal ravine fan", "CHAÑARAL" = "Salado valley, coast",
             "DIEGO DE ALMAGRO" = "Salado valley", "TIERRA AMARILLA" = "Copiapó valley",
             "COPIAPO" = "Copiapó valley", "ALTO DEL CARMEN" = "Huasco tributary fans")
t$etiqueta <- sprintf("%s — %s\n(%d / %d blocks)", nombres[t$comuna], entorno[t$comuna], t$afectadas, t$n)
t$rol <- factor(ifelse(t$comuna == "TALTAL", "Independent\nvalidation", "Calibration\ncommunes"),
                c("Calibration\ncommunes", "Independent\nvalidation"))
t$metrica <- factor(t$predictor, c("HAND", "distance to channel", "elevation", "slope"),
                    c("HAND", "Distance to channel", "Elevation", "Slope"))
orden <- unique(t$etiqueta[order(t$comuna == "TALTAL", -t$afectadas)])
t$etiqueta <- factor(t$etiqueta, rev(orden))

fig <- ggplot(t, aes(auc, etiqueta, colour = metrica, shape = metrica)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "22", colour = "grey60") +
  geom_errorbar(aes(xmin = ic95_inf, xmax = ic95_sup), width = 0.2, linewidth = 0.45,
                orientation = "y", position = position_dodge(width = 0.7)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.7), fill = "white", stroke = 0.7) +
  facet_grid(rol ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("HAND" = "#0072B2", "Distance to channel" = "#CC79A7",
                                 "Elevation" = "#E69F00", "Slope" = "#999999"), name = NULL) +
  scale_shape_manual(values = c("HAND" = 16, "Distance to channel" = 17, "Elevation" = 15, "Slope" = 18),
                     name = NULL) +
  scale_x_continuous("AUC (affected vs. non-affected census blocks)", limits = c(0.1, 1),
                     breaks = seq(0.1, 1, 0.1), expand = expansion(mult = 0.01)) +
  scale_y_discrete(NULL) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5),
        axis.text.y = element_text(size = 7, lineheight = 0.9),
        legend.position = "top", legend.box = "horizontal",
        panel.spacing = unit(6, "pt"))

save_fig(fig, "fig3_validacion_externa", width_cm = 14, height_cm = 13)
