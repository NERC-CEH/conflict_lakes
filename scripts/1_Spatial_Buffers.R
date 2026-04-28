#Script 1:
#• Pre-process datasets and locate centroids from hydrolake polys
#• Maybe also bring in catchment polys?
#• Filter UCDP conflict data by year and map both lake/UCDP datasets by year
#• Proximity analysis - generate circular buffers around conflict points (10,25,50,100km) and use st_intersect to figure out which lake centroids fall in conflict buffers
#• Analysis and plots to see if conflict exposure changes based on different buffer differences?
#• Run this analysis for particular 'anchor' years (yearly? Every 5 years?) to get time series data

# 1. Libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# 2. Load data ----
load("data/globo_topo_poly.rda") # catchment polygons
conflict <- read_csv("data/GEDEvent_v25_1.csv")

# 3. Pre-process data ----
# convert conflict to sf points
conflict_sf <- conflict %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(54009)

# calculate centroids
lake_centroids <- globo_topo_poly %>%
  st_make_valid() %>%
  st_centroid() %>%
  st_transform(54009)


buffers <- c(10000, 25000, 50000, 100000) # metres
anchor_years <- seq(2000, 2024, by = 1)   # every year

results_list <- list()

for (yr in anchor_years) {
  message(paste("Processing year:", yr))

  #  different count of conflicts per year
  conflict_yr <- conflict_sf %>% filter(year == yr)

  # seperate out by buffer differences
  year_buffer_results <- map_df(buffers, function(dist) {

    # calctualte buffers around conflict points
    conf_buffer <- st_buffer(conflict_yr, dist) %>%
      st_union()

    #  which lake centroids fall within conflict zones TRUE or FALSE
    exposed <- st_intersects(lake_centroids, conf_buffer, sparse = FALSE)

    #  summary stats distance/year
    data.frame(
      year = yr,
      buffer_km = dist / 1000,
      total_lakes = nrow(lake_centroids),
      exposed_count = sum(exposed),
      exposure_pct = (sum(exposed) / nrow(lake_centroids)) * 100
    )
  })

  results_list[[as.character(yr)]] <- year_buffer_results
}

# remove x = 2025
exposure_ts <- exposure_ts %>% filter(year < 2025)

# 4. Map and Exposure Plot ----
exposure_ts <- bind_rows(results_list)

# time series for sensitivity
exposure_buffer <- ggplot(exposure_ts, aes(x = year, y = exposure_pct, color = as.factor(buffer_km))) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Global Lake Exposure to Conflict (2000-2024)",
       x = "Year", y = "% of Global lakes exposed to conflict",
       color = "Buffer distance of lake centroid to conflict (km)") +
  theme_minimal()


# save plot
ggsave("output/exposure_sensitivity_plot.png", exposure_buffer, width = 10, height = 6)


# save exposure data
save(exposure_ts, lake_centroids, file = "data/processed_exposure_data.RData")

# get basemap
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(54009)

# 2024 and 50km
conflict_2024_buffers <- conflict_sf %>%
  filter(year == 2024) %>%
  st_buffer(50000)

# join with lake centroids to get exposed lakes in 2024
exposed_lakes_2024 <- st_join(lake_centroids, conflict_2024_buffers, left = FALSE) %>%
  mutate(violence_label = case_when(
    type_of_violence == 1 ~ "State-based",
    type_of_violence == 2 ~ "Non-state",
    type_of_violence == 3 ~ "One-sided",
    TRUE ~ "Other"
  ))

# map plot
map_2024 <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray80", size = 0.2) +
  geom_sf(data = exposed_lakes_2024,
          aes(color = violence_label, size = best),
          alpha = 0.5) +
  scale_color_manual(values = c(
    "State-based" = "orange",
    "Non-state" = "purple",
    "One-sided" = "#4DAF4A")) +
  scale_size_continuous(range = c(0.5, 4), breaks = c(10, 100, 1000)) +
  labs(title = "Lakes Exposed to Armed Conflict (2024)",
       subtitle = "Lakes within 50km of conflict events",
       color = "Type of Violence",
       size = "Fatality Estimate") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA),
    legend.position = "bottom")

# ggsave
ggsave("output/exposed_lakes_2024_map.png", map_2024, width = 12, height = 8)

# 5. Count of conflicts by country and year ----
conflict_counts <- conflict %>%
  group_by(year, country) %>%
  summarise(conflict_count = n(), .groups = "drop")
# nice table for lakes near conflicts 2024
lakes_near_conflict_2024 <- exposed_lakes_2024 %>%
  st_drop_geometry() %>%
  group_by(country) %>%
  summarise(lakes_near_conflict = n(), .groups = "drop") %>%
  arrange(desc(lakes_near_conflict)) %>%
  slice_head(n = 10)

# unique lake counts per country
lakes_near_conflict_2024 <- exposed_lakes_2024 %>%
  st_drop_geometry() %>%
  distinct(D_hylak_id, .keep_all = TRUE) %>%
  group_by(country) %>%
  summarise(unique_lakes_exposed = n(), .groups = "drop") %>%
  arrange(desc(unique_lakes_exposed)) %>%
  slice_head(n = 10)

# conflict events per country
conflict_event_counts_2024 <- conflict_sf %>%
  filter(year == 2024) %>%
  st_drop_geometry() %>%
  group_by(country) %>%
  summarise(event_count = n(), .groups = "drop") %>%
  arrange(desc(event_count))


