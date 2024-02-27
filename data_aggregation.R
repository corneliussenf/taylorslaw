
# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(sf)
library(fasterize)
library(exactextractr)

### Get data

disturbance_dat <- list.files("rawdata", pattern = glob2rx("*biotic*.csv"), full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

# Define bounding box
xmn <- min(disturbance_dat$x)
xmx <- max(disturbance_dat$x)
ymn <- min(disturbance_dat$y)
ymx <- max(disturbance_dat$y)

bbox <- st_bbox(c(xmin = xmn, xmax = xmx, ymax = ymn, ymin = ymx), crs = st_crs(3035))

# ecoregion <- read_sf("data/ecoregions.gpkg")
# ecoregion <- ecoregion %>%
#   st_transform(., crs = st_crs(3035)) %>%
#   st_crop(., bbox)

### Loop through resolutions

resolutions <- c(
  10000,
  20000,
  40000,
  80000,
  160000
)

for (res in resolutions) {
  
  print(res)
  
  ### Define grid
  
  grid <- st_make_grid(st_as_sfc(bbox), cellsize = res)
  grid <- grid %>% 
    st_as_sf() %>% 
    mutate(id = 1:n()) %>%
    rename(geom = x)
  
  write_sf(grid, paste0("data/grid_", res, ".gpkg"))
  
  ### Aggregate disturbances

  disturbance_dat_sf <- disturbance_dat %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(3035))

  disturbance_dat_grid_sf <- st_intersection(grid, disturbance_dat_sf)

  disturbance_dat_grid_sf_summary <- disturbance_dat_grid_sf %>%
    st_drop_geometry() %>%
    rename(agent_grouped = prediction_class) %>%
    group_by(id, year, agent_grouped) %>%
    summarize(disturbance_m2 = sum(area))

  ### Aggregate forest

  countries <- list.files("../attribution/results/update/predictions", pattern = glob2rx("*biotic*.csv")) %>%
    gsub("prediction_including_biotic_", "", .) %>%
    gsub(".csv", "",  .)

  grid_select <- grid %>%
    filter(id %in% unique(disturbance_dat_grid_sf_summary$id))

  fore_collector <- vector("list", length(countries))
  k <- 0

  for (c in countries) {
    print(c)
    k <- k + 1
    fore <- raster(paste0("../mapping/results/version1.0/", c, "/prediction_forestcover_", c, ".tif"))
    fore <- exact_extract(fore, grid_select, fun = "sum")
    fore_collector[[k]] <- data.frame(id = grid_select$id, forest_m2 = fore * 900) %>% filter(forest_m2 > 0)
  }

  forest <- fore_collector %>%
    bind_rows()

  forest <- forest %>%
    group_by(id) %>%
    summarize(forest_m2 = sum(forest_m2))
  
  # Temporal aggregation
  
  dat_annual <- disturbance_dat_grid_sf_summary %>%
    ungroup() %>%
    left_join(forest, by = "id") %>%
    mutate(rate = disturbance_m2 / forest_m2) %>%
    filter(!is.na(agent_grouped)) %>%
    filter(!is.na(rate)) %>%
    filter(rate <= 1) %>%
    group_by(id, agent_grouped) %>%
    summarise(mean = mean(rate),
              median = median(rate),
              var = var(rate),
              sd = sd(rate),
              iqr = IQR(rate),
              max = max(rate),
              mean_log = mean(log(rate)),
              sd_log = sd(log(rate)),
              var_log = var(log(rate)),
              median_log = median(log(rate)),
              iqr_log = IQR(log(rate)),
              mean_sqrt = mean(sqrt(rate)),
              sd_sqrt = sd(sqrt(rate)),
              var_sqrt = var(sqrt(rate)),
              median_sqrt = median(sqrt(rate)),
              iqr_sqrt = IQR(sqrt(rate)))
  
  write_csv(dat_annual, 
            paste0("data/disturbance_summary_res", res, ".csv"))
  
}

