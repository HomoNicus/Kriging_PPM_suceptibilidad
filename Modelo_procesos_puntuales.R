# ==============================================================================
# TRABAJO 2: MODELAMIENTO ESPACIAL
# ==============================================================================

# 0. LIBRERIAS Y PREPARACIÓN ----
library(tidyverse)    
library(sf)           
library(raster)       
library(gstat)        
library(spatstat)     
library(viridis)      
library(sp)           

# --- Configuración de Rutas ---
dir_base <- "C:/Users/Francisco/OneDrive - Universidad Austral de Chile/Escritorio/Estadistica_Espacial/Trabajo2"
dir_shp  <- file.path(dir_base, "Landslide_Daunia")
dir_krig <- file.path(dir_base, "resultados_kriging") 

# --- Carga de Datos ---
landslide_sf <- st_read(file.path(dir_shp, "Landslide_Daunia.shp"), quiet = TRUE)
r_elev_pred  <- raster(file.path(dir_krig, "_krigingDEM_pred.tif"))
r_slope_pred <- raster(file.path(dir_krig, "_krigingSLOPE_pred.tif"))

# --- Proyecciones ---
if (st_crs(landslide_sf)$wkt != st_crs(r_elev_pred)$wkt) {
  landslide_sf <- st_transform(landslide_sf, crs = st_crs(r_elev_pred))
}

# --- Conversión a Spatstat ---
to_im <- function(r) {
  m <- as.matrix(r); m <- m[nrow(m):1, , drop=FALSE]
  im(m, xcol=xFromCol(r, 1:ncol(r)), yrow=yFromRow(r, 1:nrow(r)))
}

Z_elev  <- to_im(r_elev_pred)
Z_slope <- to_im(r_slope_pred)
W <- as.owin(Z_elev)
coords <- st_coordinates(landslide_sf)
pp_final <- ppp(coords[,1], coords[,2], window=W, checkdup=TRUE)[W]

# --- Reescalado a Km  ---
pp_km      <- rescale(pp_final, 1000, "km")
Z_elev_km  <- rescale(Z_elev, 1000, "km")
Z_slope_km <- rescale(Z_slope, 1000, "km")
covs_km <- list(elev = Z_elev_km, slope = Z_slope_km)


# MODELAMIENTO DE LA INTENSIDAD NO-HOMOGÉNEA

# A) Modelo Poisson Base (Solo Covariables Ambientales)
fit_env   <- ppm(pp_km ~ elev + slope, covariates = covs_km)

# B) Modelo con Tendencia Lineal (Coordenadas X + Y) -> Para Figura 2 Panel B
fit_trend <- ppm(pp_km ~ elev + slope + x + y, covariates = covs_km)

# C) Modelo con Tendencia Polinómica (Curvatura) -> Para Figura 2 Panel C
fit_poly  <- ppm(pp_km ~ elev + slope + polynom(x, y, 2), covariates = covs_km)

# D) Modelo Cluster (kPPM Thomas) - Intensidad Ambiental + Interacción
fit_cluster <- kppm(pp_km ~ elev + slope, clusters = "Thomas", covariates = covs_km, method = "palm") 

# Evaluación de Significancia Estadística ----
cat("\n>>> REPORTE DE SIGNIFICANCIA (COEFICIENTES) <<<\n")

cat("\n--- A) Modelo Poisson (PPM - Solo Ambiente) ---\n")
print(summary(fit_env)$coefs.SE.CI)

cat("\n--- D) Modelo Cluster (kPPM)  ---\n")
# Cálculo para Cluster
beta_kppm <- coef(fit_cluster)
se_kppm <- sqrt(diag(vcov(fit_cluster)))
z_kppm <- beta_kppm / se_kppm
p_kppm <- 2 * pnorm(abs(z_kppm), lower.tail = FALSE)

res_kppm <- data.frame(
  Estimate = round(beta_kppm, 4), 
  Std.Error = round(se_kppm, 4), 
  Z.value = round(z_kppm, 2), 
  P.value = format.pval(p_kppm, eps=0.001), 
  Signif = ifelse(p_kppm < 0.001, "***", ifelse(p_kppm < 0.05, "*", "ns"))
)
print(res_kppm)

# Generación de Mapas Predichos de Intensidad ----
# Función auxiliar para graficar en Metros
pred_to_meters_df <- function(model) {
  pred_km <- predict(model, type="trend")
  pred_m  <- rescale(pred_km, 1/1000, "m") # Km -> Metros
  return(as.data.frame(pred_m))
}

# --- FIGURA 1: Modelo con Tendencia Polinómica ---
df_poly <- pred_to_meters_df(fit_poly)
gg_f1 <- ggplot(df_poly, aes(x, y, fill=value)) + geom_raster() +
  scale_fill_viridis(option="magma", direction=-1, name="Intensidad") + theme_light() +
  labs(title="Modelo con Tendencia Polinómica")
ggsave("Figura1_Intensidad_Polinomica.png", plot=gg_f1, width=8, height=7, dpi=300)

# --- FIGURA 2: Comparación de Ambiente vs Lineal vs Polinómica ---
df_comp <- bind_rows(
  data.frame(pred_to_meters_df(fit_env),   Modelo="A) Covariables Ambientales"),
  data.frame(pred_to_meters_df(fit_trend), Modelo="B) Tendencia Lineal (x + y)"),
  data.frame(pred_to_meters_df(fit_poly),  Modelo="C) Tendencia Polinómica")
)

gg_f2 <- ggplot(df_comp, aes(x, y, fill=value)) + 
  geom_raster() + 
  facet_wrap(~Modelo) +
  scale_fill_viridis(option="magma", direction=-1, name="Eventos/m²") + 
  theme_light() + 
  coord_fixed() +
  labs(title="Figura 2: Comparación de Estructuras de Intensidad", 
       subtitle="Evaluación de componentes de tendencia (Ambiente vs Geometría)",
       x="Este (m)", y="Norte (m)")

ggsave("Figura2_Comparacion_Modelos_3Paneles.png", plot=gg_f2, width=12, height=5, dpi=300)

# --- FIGURA 3: Mapas Finales (PPM vs kPPM) en Metros ---
# 3A: PPM
df_ppm <- pred_to_meters_df(fit_env)
gg_3a <- ggplot(df_ppm, aes(x, y, fill=value)) + geom_raster() + coord_fixed() + theme_bw() +
  scale_fill_viridis(option="magma", direction=-1, name="Eventos/m²") +
  labs(title="Intensidad Modelo Poisson (PPM)", x="Este (m)", y="Norte (m)")
ggsave("Figura3A_Mapa_PPM.png", plot=gg_3a, width=8, height=7, dpi=300)

# 3B: Cluster
df_clus <- pred_to_meters_df(fit_cluster)
gg_3b <- ggplot(df_clus, aes(x, y, fill=value)) + geom_raster() + coord_fixed() + theme_bw() +
  scale_fill_viridis(option="magma", direction=-1, name="Eventos/m²") +
  labs(title="Intensidad Modelo Cluster (kPPM)", x="Este (m)", y="Norte (m)")
ggsave("Figura3B_Mapa_Cluster.png", plot=gg_3b, width=8, height=7, dpi=300)


# ANÁLISIS DE DEPENDENCIA ESPACIAL RESIDUAL
# Calculo y Visualización de Residuos (CORREGIDO CON HARMONISE) ----
cat(">>> Calculando Residuos (Metodología: Densidad Observada - Intensidad Predicha)...\n")

# Definir el "Patrón Observado" 
densidad_obs_im <- density(pp_km, sigma=1.0)

# Función para calcular y graficar la diferencia (Residuo)
plot_residuals_methodological <- function(model, title, filename) {
  
  # Obtener la "Intensidad Predicha" por el modelo
  lambda_pred_im <- predict(model, type="trend")
  # Forzamos a que ambas imágenes tengan exactamente la misma grilla de píxeles
  imagenes_ajustadas <- harmonise(densidad_obs_im, lambda_pred_im)
  densidad_final <- imagenes_ajustadas[[1]]
  prediccion_final <- imagenes_ajustadas[[2]]
  
  # Cálculo del Residuo Espacial
  residuos_im <- eval.im(densidad_final - prediccion_final)
  

    res_m <- rescale(residuos_im, 1/1000, "m") 
  df_res <- as.data.frame(res_m)
  
  gg_res <- ggplot(df_res, aes(x=x, y=y, fill=value)) + 
    geom_raster() +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, 
                         name="Diferencia\n(Obs - Pred)") +
    coord_fixed() + theme_bw() +
    labs(title=title, 
         x="Este (m)", y="Norte (m)")
  
  ggsave(filename, plot=gg_res, width=8, height=7, dpi=300)
}

# --- Generar Mapa de Residuos para Modelo Poisson (PPM) ---
plot_residuals_methodological(fit_env, 
                              "Residuos Espaciales (Modelo Poisson)", 
                              "Figura4A_Residuos_PPM.png")

# --- Generar Mapa de Residuos para Modelo Cluster (kPPM) ---
plot_residuals_methodological(fit_cluster, 
                              "Residuos Espaciales (Modelo Cluster)", 
                              "Figura4B_Residuos_Cluster.png")

# Análisis de Agrupamiento Residual (Envelopes) ----

# --- FIGURA 5: Diagnóstico PPM (Poisson) ---
cat("   -> Generando Figura 5 (Envelope Poisson - 19 sim)...\n")
L_ppm_env <- envelope(fit_env, fun=Linhom, nsim=19, rank=1, rmax=5.0, 
                      correction="border", global=FALSE, verbose=TRUE)

png("Figura5_Envelope_PPM.png", width=2000, height=2000, res=300)
par(mar=c(5,5,4,2))
plot(L_ppm_env, . - r ~ r, main="Dependencia Residual (PPM)", 
     legend=TRUE, shade=c("hi", "lo"), ylab="L_inhom(r) - r")
abline(h=0, lty=2, col="red")
dev.off()

# --- FIGURA 6: Diagnóstico kPPM (Cluster) ---
cat("   -> Generando Figura 6 (Envelope Cluster - 19 sim)...\n")
L_clus_env <- envelope(fit_cluster, fun=Linhom, nsim=19, rank=1, rmax=5.0, 
                       correction="border", global=FALSE, verbose=TRUE)

png("Figura6_Envelope_Cluster.png", width=2000, height=2000, res=300)
par(mar=c(5,5,4,2))
plot(L_clus_env, . - r ~ r, main="Dependencia Residual (Cluster)", 
     legend=TRUE, shade=c("hi", "lo"), ylab="L_inhom(r) - r")
abline(h=0, lty=2, col="gray50")
dev.off()

# COMPARACIÓN VISUAL DE SIMULACIONES
png("Figura7_Comparacion_Simulaciones.png", width=3000, height=1000, res=300)
par(mfrow=c(1,3), mar=c(1,1,3,1))

# A) Datos Reales
plot(pp_km, pch=16, cex=0.5, cols="black", main="A) Datos Reales", use.marks=F)

# B) Simulación Poisson (Dispersa)
plot(simulate(fit_env, nsim=1, drop=TRUE), pch=16, cex=0.5, cols="blue", 
     main="B) Simulación Poisson", use.marks=F)

# C) Simulación Cluster (Agrupada)
plot(simulate(fit_cluster, nsim=1, drop=TRUE), pch=16, cex=0.5, cols="red", 
     main="C) Simulación Cluster", use.marks=F)

dev.off()

# Comparación AIC
cat(paste("AIC Poisson:", round(AIC(fit_env), 2), "\n"))
cat(paste("AIC Cluster:", round(AIC(fit_cluster), 2), "\n"))

# Parámetros del Cluster
pars <- fit_cluster$par
if(is.null(pars)) pars <- fit_cluster$modelpar
cat("\nParámetros de Estructura Espacial (Cluster Thomas):\n")
print(round(pars, 4))
cat("(Kappa = Densidad de clusters, Sigma = Radio de dispersión)\n")
