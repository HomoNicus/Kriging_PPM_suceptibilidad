library(tidyverse)
library(sf)    
library(raster)
library(units)


landslide_sf <- st_read(
  "F:/Estadistica_DataSet/GeomorphologicalLandslide/centroides/Landslide_Daunia.shp"
)

dem <- raster("F:/Estadistica_DataSet/SRTM_30m/Altitud_AE.tif")
slope <- raster("F:/Estadistica_DataSet/SRTM_30m/Pendiente_AE.tif")


# transformamos a dataframe los landslide y los raster ----
# deslizamientos
landslide_df <- landslide_sf %>%
  st_drop_geometry() %>%
  as_tibble()

glimpse(landslide_df)

#Altitud

dem_df <- as.data.frame(dem, xy = TRUE) %>%
  as_tibble() %>%
  rename(elevation = 3)

summary(dem_df$elevation)
#Pendiente
slope_df <- as.data.frame(slope, xy = TRUE) %>%
  as_tibble() %>%
  rename(slope = 3)

# Agregamos coordenadas al dataframe (a nuestros deslizamientos)
coords <- st_coordinates(landslide_sf)

landslide_df_xy <- landslide_sf %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  )


# Visualizamos nuestos datos----
glimpse(landslide_df_xy)


# Ploteamos nuestros datos ----
library(ggplot2)

ggplot(landslide_df_xy, aes(x = x, y = y)) +
  geom_point(aes(color = Pendiente_), size = 1.4, alpha = 0.8) +
  coord_equal() +
  scale_color_viridis_c(
    name = "Pendiente (°)",
    option = "C"
  ) +
  labs(
    title = "Deslizamientos en Daunia",
    subtitle = "Coloreados según pendiente",
    x = "Este (m)",
    y = "Norte (m)"
  ) +
  theme_minimal()



# Visualizamos distribucion de puntos y distribucion espacial de variables continuas----
# Es decir, visualizamos en dos mapas la pendiente y la altitud junto con los puntos (deslizamientos)
library(patchwork)
#Agrego la geometria a los puntos
landslide_Daunia_crs <- st_as_sf(landslide_df_xy, coords = c("x", "y"), crs = 32633)

# 1. Mapa de altitud con puntos de deslizamientos
p1 <- ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_viridis_c(name = "Altitud (m)") +
  geom_sf(data = landslide_Daunia_crs, color = "red", size = 0.25, alpha = 0.7) +
  coord_sf(crs = st_crs(landslide_Daunia_crs)) +
  theme_minimal() +
  labs(title = "Deslizamientos sobre mapa de altitud",
       x = "Longuitud",
       y = "Latitud")

# 2. Mapa de pendiente con puntos de deslizamientos
p2 <- ggplot() +
  geom_raster(data = slope_df, aes(x = x, y = y, fill = slope)) +
  scale_fill_viridis_c(name = "Pendiente (°)") +
  geom_sf(data = landslide_Daunia_crs, color = "red", size = 0.25, alpha = 0.7) +
  coord_sf(crs = st_crs(landslide_Daunia_crs)) +
  theme_minimal() +
  labs(title = "Deslizamientos sobre mapa de pendiente",
       x = "Longuitud",
       y = "Latitud")

# 3. Mostrar en paralelo
p1 + p2  # patchwork los muestra uno al lado del otro



