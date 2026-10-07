# Fig. 2 — Validation against the March 2015 debris flow in Taltal (unit level).
# (a) ROC curves; (b) AUC (95 % bootstrap CI) and Spearman rho with affected population.
# Ejecutar desde la raíz del proyecto: Rscript paper/R/fig2_validacion_taltal.R
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(pROC) })
source("paper/R/theme_paper.R")
setup_paper_theme("elsevier")

d <- targets::tar_read(tabla_2015, store = st_corr)
d$riesgo_final <- d$p_alcance * (d$v_extendida_indice + d$exposicion + d$locc) / 3
d$riesgo_aditivo <- 0.25 * (d$susceptibilidad + d$v_extendida_indice + d$exposicion + d$locc)
modelos <- c("Initiation susceptibility" = "susceptibilidad", "Conventional additive index" = "riesgo_aditivo",
             "Reach index" = "p_alcance", "Reach-based risk (this study)" = "riesgo_final")

set.seed(42)
curvas <- list(); resumen <- list()
for (n in names(modelos)) {
  r <- roc(d$afectada, d[[modelos[[n]]]], direction = "<", levels = c(FALSE, TRUE), quiet = TRUE)
  ci <- ci.auc(r, method = "bootstrap", boot.n = 2000)
  rho <- cor(d[[modelos[[n]]]], d$pob_afectada, method = "spearman")
  curvas[[n]] <- data.frame(modelo = n, fpr = 1 - r$specificities, tpr = r$sensitivities)
  resumen[[n]] <- data.frame(modelo = n, auc = ci[2], lo = ci[1], hi = ci[3], rho = rho)
}
curvas <- do.call(rbind, curvas); resumen <- do.call(rbind, resumen)
curvas$modelo <- factor(curvas$modelo, names(modelos)); resumen$modelo <- factor(resumen$modelo, names(modelos))
etiqueta <- setNames(sprintf("%s  (AUC %.2f)", names(modelos), resumen$auc), names(modelos))

pa <- ggplot(curvas[order(curvas$modelo, curvas$fpr, curvas$tpr), ], aes(fpr, tpr, colour = modelo)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "22", colour = "grey60") +
  geom_step(linewidth = 0.6, direction = "hv") +
  scale_colour_manual(values = col_modelo, name = NULL) +
  scale_x_continuous("False positive rate", limits = c(0, 1), expand = expansion(mult = 0.01)) +
  scale_y_continuous("True positive rate", limits = c(0, 1), expand = expansion(mult = 0.01)) +
  coord_equal() +
  guides(colour = guide_legend(nrow = 1))

largo <- rbind(
  data.frame(modelo = resumen$modelo, metrica = "AUC (affected / not affected)", valor = resumen$auc,
             lo = resumen$lo, hi = resumen$hi),
  data.frame(modelo = resumen$modelo, metrica = "Spearman \u03c1 (affected population)", valor = resumen$rho,
             lo = NA, hi = NA))
ref <- data.frame(metrica = unique(largo$metrica), x = c(0.5, 0))
pb <- ggplot(largo, aes(valor, modelo, colour = modelo)) +
  geom_vline(data = ref, aes(xintercept = x), linewidth = 0.3, linetype = "22", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.5, orientation = "y", na.rm = TRUE) +
  geom_point(size = 2) +
  facet_wrap(~metrica, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = col_modelo, guide = "none") +
  scale_y_discrete(NULL, limits = rev(names(modelos)),
                   labels = function(x) gsub(" \\(", "\n(", x)) +
  scale_x_continuous(NULL) +
  theme(strip.text = element_text(hjust = 0, size = 8), axis.text.y = element_text(size = 7),
        panel.spacing = unit(8, "pt"))

fig <- (pa | pb) + plot_layout(widths = c(1, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom", legend.key.width = unit(14, "pt"),
        legend.text = element_text(size = 7), legend.margin = margin(0, 0, 0, 0))
save_fig(fig, "fig2_validacion_taltal", width_cm = 18, height_cm = 10)
print(resumen, row.names = FALSE)
