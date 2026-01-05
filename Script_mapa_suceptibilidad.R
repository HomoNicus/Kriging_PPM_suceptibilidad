
library(raster)
library(ggplot2)

# Ruta al raster de intensidad PPM ----
path_ppm <- "~/Ciencia de datos/Espacial/Trabajo Final/Mapa_Intensidad_Predicha_Cluster/Mapa_Intensidad_Predicha_PPM.tif"

#  Función: clasificación por cuantiles ----
quantile_classify <- function(r, probs = seq(0, 1, by = 0.2)) {
  v_all <- getValues(r)
  v_ok  <- v_all[is.finite(v_all)]
  
  brks <- quantile(v_ok, probs = probs, na.rm = TRUE)
  brks <- sort(unique(brks))
  
  k <- length(brks) - 1
  if (k < 2) stop("No hay variabilidad suficiente para clasificar.")
  
  cl <- findInterval(
    v_all,
    vec = brks,
    rightmost.closed = TRUE,
    all.inside = TRUE
  )
  cl[!is.finite(v_all)] <- NA
  
  r_class <- setValues(raster(r), cl)
  
  list(r_class = r_class, breaks = brks, k = k)
}

#  Función de ploteo ----
plot_susc <- function(r_class, title, k) {
  cols <- c("green4", "green2", "yellow", "orange", "red3")[seq_len(k)]
  labs_k <- c("Muy baja", "Baja", "Media", "Alta", "Muy alta")[seq_len(k)]
  
  df <- as.data.frame(r_class, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "clase"
  df$clase <- factor(df$clase, levels = seq_len(k), labels = labs_k)
  
  ggplot(df, aes(x = x, y = y, fill = clase)) +
    geom_raster() +
    coord_equal() +
    scale_fill_manual(values = cols, drop = FALSE) +
    theme_light() +
    labs(
      title = title,
      fill  = "Susceptibilidad",
      x     = "Este (m)",
      y     = "Norte (m)"
    )
}

#  Leer raster PPM ----
r_ppm <- raster(path_ppm)

# Clasificar ----
q_ppm <- quantile_classify(r_ppm)

#  Plot final ----
p_ppm <- plot_susc(
  q_ppm$r_class,
  "Mapa de Susceptibilidad Espacial (PPM)",
  q_ppm$k
)

print(p_ppm)