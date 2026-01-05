# Cargar librerías necesarias
library(sf)          # manejo de datos espaciales (shapefiles)
library(tidyverse)   # incluye purrr, dplyr, ggplot2, etc.
library(rnaturalearth)      # datos de Natural Earth
library(rnaturalearthdata)  # datasets de Natural Earth
library(ggplot2)

# 1. Listar y leer los shapefiles que comienzan con "lim"
shp_files <- list.files(path = "F:/Estadistica_DataSet/GeomorphologicalLandslide/GeomorphologicalLandslide/GeomorphologicalLandslideInventoryMap/GeomorphologicalLandslideInventoryMap", pattern = "^lim.*\\.shp$", full.names = TRUE)
shp_list  <- map(shp_files, st_read)  # leer cada shapefile como sf

# 2. Asegurarse de que todos están en CRS EPSG:32633
shp_list <- map(shp_list, ~ st_transform(.x, 32633))

# 3. Calcular centroids de cada polígono
centroid_list <- map(shp_list, st_centroid)



# 4. Unificar columnas: obtener la lista completa de nombres de columnas
all_cols <- reduce(map(centroid_list, names), union)

# Agregar columnas faltantes en cada sf y reordenar columnas según all_cols
centroid_list <- map(centroid_list, function(sf_obj) {
  missing_cols <- setdiff(all_cols, names(sf_obj))
  if(length(missing_cols) > 0) {
    sf_obj[missing_cols] <- NA   # añade columnas faltantes con NA
  }
  sf_obj <- sf_obj[, all_cols]   # reordena las columnas según all_cols
  return(sf_obj)
})


# 5. Combinar (apilar) todos los centroides en un solo sf
all_centroids <- do.call(rbind, centroid_list)

# Verificar CRS y columnas comunes
st_crs(all_centroids)$epsg   # debería ser 32633
names(all_centroids)        # columnas disponibles en el objeto combinado


# 6. Cargar un mapa base mundial y transformarlo a EPSG:32633
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, 32633)

# 7. Graficar con ggplot2
# Calcular límites (bounding box) de los centroides
bbox <- st_bbox(all_centroids)

# Visualizar el mapa enfocado automáticamente
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "white") +
  geom_sf(data = all_centroids, aes(color = "red"), size = 2, alpha = 0.7) +
  coord_sf(
    crs = st_crs(all_centroids),
    xlim = c(bbox$xmin, bbox$xmax),
    ylim = c(bbox$ymin, bbox$ymax),
    expand = FALSE
  ) +
  theme_light() +
  labs(
    title = "Centroides de deslizamientos - Región de Daunia (Italia)",
    color = "Generación"
  )

# Exportar los centroides como shapefile
st_write(all_centroids, "F:/Estadistica_DataSet/GeomorphologicalLandslide/centroides/centroides_daunia.shp", delete_dsn = TRUE)

