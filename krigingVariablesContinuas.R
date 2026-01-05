
# Metodología paso a paso para realizar kriging ----

library(raster)
library(sp)
library(gstat)
library(tidyverse)

set.seed(123)

# 1) Cargar y preprocesar datos de pendiente y elevacion -----

dem   <- raster("F:/Estadistica_DataSet/SRTM_30m/Altitud_AE.tif")
slope <- raster("F:/Estadistica_DataSet/SRTM_30m/Pendiente_AE.tif")

# Elegir variable a analizar:
r <- dem
varname <- "DEM"   # cambia a "SLOPE" y r <- slope si quieres

# Inspección básica de las variables a interpolar
print(r)
cat("\nCRS:\n"); print(crs(r))
cat("\nResolución:\n"); print(res(r))
cat("\nExtent:\n"); print(extent(r))

# Visualizacion rapida de raster cargado
plot(r, main = paste0(varname, " (raster original)"))

# Parámetros principales: Ajuste de parametros principales del tipo de grilla a interpolar y los
# modelos a evaluar (esferico, exponencial y gaussiano)
n_sample <- 3000
use_trend <- TRUE               # TRUE -> Universal: z ~ x + y; FALSE -> Ordinario: z ~ 1
grid_fact <- 3                  
models_to_try <- c("Sph","Exp","Gau")

# Opción anisotropía:
use_anisotropy <- TRUE          # TRUE = calcula variogramas direccionales y permite ajuste anisotrópico
alphas <- c(0,45,90,135)
tol.hor <- 22.5

# 2) Puntos de muestreo (extraer observaciones desde el raster) ----
# Realizamos malla de puntos de muestreo aleatoreo en nuestro raster de muestra. 
samp_df <- raster::sampleRandom(r, size = n_sample, xy = TRUE, na.rm = TRUE) %>% as.data.frame()
names(samp_df) <- c("x","y","z")

# TRealizamos tabla resumida para el informe
summary_tbl <- samp_df %>%
  summarise(
    n = n(),
    z_min = min(z, na.rm = TRUE),
    z_mean = mean(z, na.rm = TRUE),
    z_max = max(z, na.rm = TRUE)
  )
print(summary_tbl)

# Plot de puntos muestreados (coloreados por z)
ggplot(samp_df, aes(x = x, y = y, color = z)) +
  geom_point(size = 0.7, alpha = 0.6) +
  coord_equal() +
  theme_minimal() +
  labs(title = paste0(varname, " - Puntos de muestreo aleatorio (n=", n_sample, ")"),
       color = varname) %>%
  print()

# Convertir a SpatialPointsDataFrame
samp_sp <- samp_df
coordinates(samp_sp) <- ~ x + y
proj4string(samp_sp) <- crs(r)

# Preparar x,y en data (para Universal Kriging)
samp_sp@data$x <- coordinates(samp_sp)[,1]
samp_sp@data$y <- coordinates(samp_sp)[,2]


# 5) Definimos la formula del kriging, en este caso es la formula del kriging universal-----

form <- if (use_trend) z ~ x + y else z ~ 1
cat("\nTipo de kriging:", if (use_trend) "Universal (z ~ x + y)" else "Ordinary (z ~ 1)", "\n")


# 3) Variograma empírico (Isotropico y anisotropico


coords <- coordinates(samp_sp)
max_dist <- max(dist(coords))
cutoff <- max_dist / 2
width  <- cutoff / 15

cat("\nParámetros variograma:\n")
cat("cutoff =", cutoff, "\n")
cat("width  =", width, "\n")

# 3A) Variograma omnidireccional (isotrópico)
v_emp <- gstat::variogram(form, samp_sp, cutoff = cutoff, width = width)

# Plot variograma empírico (ggplot)
ggplot(v_emp, aes(x = dist, y = gamma)) +
  geom_point() + geom_line() +
  theme_minimal() +
  labs(title = paste0(varname, " - Variograma empírico omnidireccional"),
       x = "Distancia", y = "Semivarianza (gamma)") %>%
  print()

# 3B) Variogramas direccionales (para diagnosticar anisotropía)
if (use_anisotropy) {
  v_dir <- lapply(alphas, function(a) {
    gstat::variogram(form, samp_sp, cutoff = cutoff, width = width, alpha = a, tol.hor = tol.hor) %>%
      mutate(alpha = a)
  }) %>% bind_rows()
  
  ggplot(v_dir, aes(x = dist, y = gamma, color = factor(alpha))) +
    geom_point(alpha = 0.6) + geom_line(alpha = 0.6) +
    theme_minimal() +
    labs(title = paste0(varname, " - Variogramas direccionales (diagnóstico anisotropía)"),
         x = "Distancia", y = "Semivarianza", color = "alpha (°)") %>%
    print()
}


# 4) Ajuste de modelo continuo de variograma (comparar Sph/Exp/Gau) -----


fit_one <- function(v_emp, model) {
  psill0 <- var(v_emp$gamma, na.rm = TRUE)
  range0 <- max(v_emp$dist, na.rm = TRUE) / 3
  nug0   <- min(v_emp$gamma, na.rm = TRUE) * 0.5
  
  vgm0 <- gstat::vgm(psill = psill0, model = model, range = range0, nugget = max(0, nug0))
  fit  <- try(gstat::fit.variogram(v_emp, vgm0), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  tibble(model = model, SSErr = attr(fit, "SSErr"), fit = list(fit))
}

fits_tbl <- map_dfr(models_to_try, ~fit_one(v_emp, .x)) %>% arrange(SSErr)

if (nrow(fits_tbl) == 0) stop("No se pudo ajustar ningún modelo.")

# Tabla para informe: ranking por SSErr
print(dplyr::select(fits_tbl, model, SSErr))  ### 

best_model <- fits_tbl$model[1]
v_fit_iso  <- fits_tbl$fit[[1]]

cat("\nMejor modelo isotrópico:", best_model, "\n")
print(v_fit_iso)

# Plot empírico + curvas de todos los modelos (comparación visual)
# (usamos variogramLine para generar líneas)
make_line_df <- function(vgm_fit, maxdist) {
  as.data.frame(gstat::variogramLine(vgm_fit, maxdist = maxdist))
}

line_all <- fits_tbl %>%
  dplyr::mutate(
    line = purrr::map(
      fit,
      ~ make_line_df(.x, maxdist = max(v_emp$dist, na.rm = TRUE))
    )
  ) %>%
  dplyr::select(model, line) %>%
  tidyr::unnest(line)

ggplot() +
  geom_point(data = v_emp, aes(x = dist, y = gamma), size = 2) +
  geom_line(data = line_all, aes(x = dist, y = gamma, color = model), linewidth = 1) +
  theme_minimal() +
  labs(title = paste0(varname, " - Variograma empírico + modelos ajustados"),
       x = "Distancia", y = "Semivarianza", color = "Modelo") %>%
  print()


# 4B) Ajuste de variograma anisotropico (solo si use_anisotropy=TRUE) ----
# aca se realiza en diagnostico de la anisotropia y la modelación de la anisotropia utilizando
# la dirección que mejor se ajuste al variograma. 


v_fit_final <- v_fit_iso  

if (use_anisotropy) {
  # Acá se ajusta el angulo y el ratio segun lo observado en el variograma direccional 
  anis_angle <- 135      # dirección principal (ej: 0,45,90,135)
  anis_ratio <- 0.5   
  
  # Crear un modelo inicial anisotrópico (misma familia que el mejor isotrópico)
  psill0 <- var(v_emp$gamma, na.rm = TRUE)
  range0 <- max(v_emp$dist, na.rm = TRUE) / 3
  nug0   <- min(v_emp$gamma, na.rm = TRUE) * 0.5
  
  vgm0_anis <- gstat::vgm(psill = psill0, model = best_model, range = range0,
                          nugget = max(0, nug0), anis = c(anis_angle, anis_ratio))
  
  v_fit_anis <- try(gstat::fit.variogram(v_emp, vgm0_anis), silent = TRUE)
  
  if (!inherits(v_fit_anis, "try-error")) {
    cat("\nModelo anisotrópico ajustado (candidato):\n")
    print(v_fit_anis)
    
    # Plot comparación isotropico vs anisotropico
    line_iso  <- make_line_df(v_fit_iso,  maxdist = max(v_emp$dist, na.rm = TRUE)) %>% mutate(type="Isotrópico")
    line_anis <- make_line_df(v_fit_anis, maxdist = max(v_emp$dist, na.rm = TRUE)) %>% mutate(type="Anisotrópico")
    
    ggplot() +
      geom_point(data = v_emp, aes(x = dist, y = gamma), size = 2) +
      geom_line(data = bind_rows(line_iso, line_anis),
                aes(x = dist, y = gamma, color = type), linewidth = 1) +
      theme_minimal() +
      labs(title = paste0(varname, " - Comparación modelo isotrópico vs anisotrópico"),
           x = "Distancia", y = "Semivarianza", color = "Modelo") %>%
      print()
    
    # Aca decidimos si utilizar el modelo anisotropico o no, segun el que mejor se ajuste. 
    v_fit_final <- v_fit_anis
  } else {
    message("No se pudo ajustar el modelo anisotrópico. Se mantiene isotrópico.")
  }
}


# 6) Validación cruzada (comparar si quieres iso vs anis) -----


run_cv_metrics <- function(v_fit) {
  cv <- gstat::krige.cv(formula = form, locations = samp_sp, model = v_fit, nfold = 10)
  df <- as.data.frame(cv)
  
  tibble(
    ME = mean(df$residual, na.rm = TRUE),
    RMSE = sqrt(mean(df$residual^2, na.rm = TRUE)),
    MSE_std = mean(df$residual / sqrt(df$var1.var), na.rm = TRUE),
    RMSSE = sqrt(mean((df$residual^2) / df$var1.var, na.rm = TRUE))
  )
}

metrics_final <- run_cv_metrics(v_fit_final)
cat("\nMétricas CV (modelo final):\n")
print(metrics_final)

# (Opcional) gráficos de diagnóstico Validacion cruzada
cv_final <- gstat::krige.cv(formula = form, locations = samp_sp, model = v_fit_final, nfold = 10)
cv_df <- as.data.frame(cv_final)

ggplot(cv_df, aes(x = observed, y = var1.pred)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0) +
  theme_light() +
  labs(title = paste0(varname, " - CV: Observado vs Predicho"),
       x = "Observado", y = "Predicho") %>%
  print()

ggplot(cv_df, aes(x = residual)) +
  geom_histogram(bins = 40) +
  theme_light() +
  labs(title = paste0(varname, " - CV: Histograma de residuos"),
       x = "Residual", y = "Frecuencia") %>%
  print()

ggplot(cv_df, aes(sample = residual)) +
  stat_qq() + stat_qq_line() +
  theme_light() +
  labs(title = paste0(varname, " - CV: QQplot de residuos")) %>%
  print()


# 7) Interpolación (kriging) + mapas predicción e incertidumbre ----


r_grid <- if (grid_fact > 1) raster::aggregate(r, fact = grid_fact) else r
grid_sp <- raster::rasterToPoints(r_grid, spatial = TRUE)
grid_sp <- as(grid_sp, "SpatialPixelsDataFrame")

system.time({
k <- gstat::krige(formula = form, locations = samp_sp, newdata = grid_sp, model = v_fit_final)
})

system.time({
r_pred <- raster(k["var1.pred"])
})

r_var  <- raster(k["var1.var"])

# Ayuda para realizar el ggplot -----
r_to_df <- function(rr, name = "value") {
  as.data.frame(rasterToPoints(rr)) %>% setNames(c("x","y",name))
}

ggplot(r_to_df(r_pred, "pred"), aes(x, y, fill = pred)) +
  geom_raster() + coord_equal() + theme_light() +
  labs(title = paste0(varname, " - Mapa Kriging (predicción)"),
       fill = "Pred") %>%
  print()

ggplot(r_to_df(r_var, "var"), aes(x, y, fill = var)) +
  geom_raster() + coord_equal() + theme_light() +
  labs(title = paste0(varname, " - Varianza de predicción (incertidumbre)"),
       fill = "Var") %>%
  print()

#  Guardamos resultados resultados -----
writeRaster(r_pred, paste0("F:/Estadistica_DataSet/Productos/_kriging", varname, "_pred.tif"), overwrite = TRUE)
writeRaster(r_var,  paste0("F:/Estadistica_DataSet/Productos/kriging_", varname, "_var.tif"),  overwrite = TRUE)
