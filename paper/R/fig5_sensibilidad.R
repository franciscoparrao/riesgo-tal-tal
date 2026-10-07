# Fig. 5 — Sensitivity of the main result.
# (a) channel-initiation threshold for HAND (recalibrated on Atacama each time) and DEM choice;
# (b) Monte Carlo weights of the consequence components (V, E, LOCC) ~ Dirichlet(1,1,1);
# (c) alternative risk formulations (AUC, 95 % CI) against the 2015 Taltal footprint.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("paper/R/theme_paper.R")
source("R/validacion_2015.R")      # auc()
setup_paper_theme("elsevier")
S <- "salidas/corregido/sensibilidad/"

# (a) umbral + DEM
u <- read.csv(paste0(S, "umbral_cauce.csv"))
dem <- read.csv(paste0(S, "dem.csv"))
ua <- rbind(data.frame(umbral = u$umbral_area_km2 <- u$area_umbral_km2, auc = u$auc_loco_media,
                       lo = NA, hi = NA, serie = "Atacama, leave-one-commune-out (mean)"),
            data.frame(umbral = u$area_umbral_km2, auc = u$auc_taltal, lo = u$ic95_inf, hi = u$ic95_sup,
                       serie = "Taltal, independent (Copernicus DEM)"))
fab <- data.frame(umbral = 1.8, auc = dem$auc[2], lo = dem$ic95_inf[2], hi = dem$ic95_sup[2],
                  serie = "Taltal, independent (FABDEM)")
cols_a <- c("Atacama, leave-one-commune-out (mean)" = "#56B4E9", "Taltal, independent (Copernicus DEM)" = "#0072B2",
            "Taltal, independent (FABDEM)" = "#009E73")
pa <- ggplot(ua, aes(umbral, auc, colour = serie)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "22", colour = "grey60") +
  annotate("rect", xmin = 1.8 * 0.93, xmax = 1.8 * 1.07, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.08, linewidth = 0.4, na.rm = TRUE) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.8) +
  geom_errorbar(data = fab, aes(ymin = lo, ymax = hi), width = 0.08, linewidth = 0.4,
                position = position_nudge(x = 0.12)) +
  geom_point(data = fab, shape = 17, size = 2, position = position_nudge(x = 0.12)) +
  scale_colour_manual(values = cols_a, name = NULL) +
  scale_x_continuous(expression("Channel initiation area (km"^2*")"), breaks = c(0.9, 1.8, 3.6)) +
  scale_y_continuous("AUC", limits = c(0.45, 0.9)) +
  guides(colour = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5), legend.key.height = unit(7, "pt"))

# (b) Monte Carlo de pesos (recalculado para tener la distribución completa)
d <- targets::tar_read(tabla_2015, store = st_corr)
set.seed(42); n <- 2000
w <- matrix(rexp(3 * n), ncol = 3); w <- w / rowSums(w)
mc <- do.call(rbind, lapply(seq_len(n), function(i) {
  s <- d$p_alcance * (w[i, 1] * d$v_extendida_indice + w[i, 2] * d$exposicion + w[i, 3] * d$locc)
  data.frame(metrica = c("AUC", "Spearman ρ"), valor = c(auc(s, d$afectada),
             cor(s, d$pob_afectada, method = "spearman")))
}))
ref_b <- data.frame(metrica = c("AUC", "Spearman ρ"),
                    final = c(auc(d$p_alcance * (d$v_extendida_indice + d$exposicion + d$locc) / 3, d$afectada),
                              cor(d$p_alcance * (d$v_extendida_indice + d$exposicion + d$locc) / 3,
                                  d$pob_afectada, method = "spearman")))
pb <- ggplot(mc, aes(valor)) +
  geom_histogram(bins = 40, fill = "#D55E00", colour = NA, alpha = 0.75) +
  geom_vline(data = ref_b, aes(xintercept = final), linewidth = 0.4, linetype = "22") +
  facet_wrap(~metrica, ncol = 1, scales = "free") +
  scale_x_continuous("AUC or Spearman ρ over 2000 random weightings") +
  scale_y_continuous("Count", expand = expansion(mult = c(0, 0.05))) +
  theme(strip.text = element_text(hjust = 0.5, size = 7.5, margin = margin(2, 0, 2, 0)))

# (c) formulaciones
f <- read.csv(paste0(S, "formulaciones_riesgo.csv"), check.names = FALSE)
etiq <- c("solo amenaza P(alcance)" = "Reach index only",
          "H x mean(V,E,L) [final, V = índice extendida]" = "Reach × mean(V,E,L)  [final]",
          "H x mean(V,E,L), V = índice base" = "Reach × mean(V,E,L), V base index",
          "H x mean(V,E,L), V = fuzzy extendida" = "Reach × mean(V,E,L), V fuzzy ext.",
          "H x mean(V,E,L), V = fuzzy base" = "Reach × mean(V,E,L), V fuzzy base",
          "H x V (extendida)" = "Reach × V",
          "H x mean(V,E) (sin LOCC)" = "Reach × mean(V,E)",
          "H x (V x E x L)^(1/3)" = "Reach × geometric mean(V,E,L)",
          "aditivo 0.25 (H,V,E,L)" = "Additive (reach,V,E,L)",
          "aditivo convencional (S,V,E,L)" = "Additive (S,V,E,L)  [conventional]")
f$etiqueta <- factor(etiq[f$formulacion], rev(etiq))
f$tipo <- ifelse(grepl("Additive", f$etiqueta), "Additive", ifelse(grepl("geometric", f$etiqueta), "Multiplicative (geometric)", "Multiplicative"))
pc <- ggplot(f, aes(auc, etiqueta, colour = tipo)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "22", colour = "grey60") +
  geom_errorbar(aes(xmin = ic95_inf, xmax = ic95_sup), width = 0.3, linewidth = 0.4, orientation = "y") +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Multiplicative" = "#D55E00", "Additive" = "#E69F00",
                                 "Multiplicative (geometric)" = "#999999"), name = NULL) +
  scale_x_continuous("AUC (95 % CI)", limits = c(0.3, 0.9)) + scale_y_discrete(NULL) +
  theme(axis.text.y = element_text(size = 7), legend.position = "none")

fig <- (free(pa | pb) / pc) + plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
save_fig(fig, "fig5_sensibilidad", width_cm = 18, height_cm = 15)
