#Script 3:
#  • Extract species richness and LPI data for each of the exposed lakes per year (only LPI can give us time series data of the two….)
#  • Plot: x axis = % lakes near conflict, y axis = biodiversity index relative to baseline
#  • Regression of conflict exposure and land use type against biodiversity metrics to see if conflict is a significant 'risk multiplier' - what about lakes DIRECTLY in a conflict area vs lakes in a country at war that is not directly in a war zone?
#  • Map hotspots where high conflict exposure overlaps with severe biodiversity decline

# 1. Load Libraries ----

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(broom)

# 2. Load Data ----
load("data/globo_topo_poly.rda")
lpd <- read_csv("data/LPD_Public_2024.csv")
sr_raster <- rast("data/IUCN_Freshwater_SR_2024.tif")

world_tags <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(name_long, region_un, continent) %>%
  st_transform(54009)
lake_centroids_tagged <- globo_topo_poly %>%
  st_make_valid() %>%
  st_centroid() %>%
  st_transform(54009) %>%
  st_join(world_tags)
lake_centroids_tagged$species_richness <- terra::extract(sr_raster, st_transform(lake_centroids_tagged, st_crs(sr_raster)))[,2]

#3.  Species Richness icun  ----

conflict_2024 <- conflict_sf %>% filter(year == 2024)
buf_10 <- st_union(st_buffer(conflict_2024, 10000))
buf_25 <- st_union(st_buffer(conflict_2024, 25000))
buf_50 <- st_union(st_buffer(conflict_2024, 50000))

war_countries_2024 <- conflict %>%
  filter(year == 2024) %>% group_by(country) %>%
  summarise(total_deaths = sum(best)) %>% filter(total_deaths >= 25) %>% pull(country)

rainfall_df <- lake_centroids_tagged %>%
  mutate(
    d10 = as.logical(st_intersects(., buf_10, sparse = FALSE)),
    d25 = as.logical(st_intersects(., buf_25, sparse = FALSE)),
    d50 = as.logical(st_intersects(., buf_50, sparse = FALSE))
  ) %>%
  st_drop_geometry() %>%
  mutate(exposure_cat = case_when(
    d10 ~ "10km Buffer",
    d25 ~ "25km Buffer",
    d50 ~ "50km Buffer",
    country %in% war_countries_2024 ~ "Country at War",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(exposure_cat), !is.na(species_richness))

rainfall_df$exposure_cat <- factor(rainfall_df$exposure_cat,
                                   levels = c("Country at War", "50km Buffer", "25km Buffer", "10km Buffer"))

ggplot(rainfall_df, aes(x = exposure_cat, y = species_richness, fill = exposure_cat)) +
  stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) +
  geom_point(size = 0.5, alpha = .1, position = position_jitter(width = .15)) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("10km Buffer"="#67000d", "25km Buffer"="#a50f15", "50km Buffer"="#ef3b2c", "Country at War"="#fc8d59")) +
  labs(title = "Species Richness in Conflict Zones (2024)", x = "Intensity", y = "Richness") +
  theme_minimal() + theme(legend.position = "none")

# 4.  Ukraine ----
ukraine_data <- lake_exposure_long %>%
  filter(year == 2024, country == "Ukraine") %>%
  inner_join(lake_metadata, by = "Hylak_id") %>%
  mutate(Group = ifelse(is_exposed, "Lakes in Conflict Zones (50km)", "All National Lakes"))

ggplot(ukraine_data, aes(x = species_richness, fill = Group)) +
  geom_density(alpha = 0.6, color = "white") +
  geom_rug(aes(color = Group), alpha = 0.3) +
  scale_fill_manual(values = c("All National Lakes" = "#0057b7",
                               "Lakes in Conflict Zones (50km)" = "#ffd700")) +
  scale_color_manual(values = c("All National Lakes" = "#0057b7",
                                "Lakes in Conflict Zones (50km)" = "#ffd700")) +
  labs( title = "Ukraine",
        x = "Raw Species Richness (IUCN 2024)",
       y = "Density (Frequency of Lakes)",
       fill = "Lake Group") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

# ggsave
ggsave("output/ukraine_species_richness_density.png", width = 10, height = 6)

# 5. Colombia ----
colombia_data <- lake_exposure_long %>%
  filter(year == 2024, country == "Colombia") %>%
  inner_join(lake_metadata, by = "Hylak_id") %>%
  mutate(Group = ifelse(is_exposed, "Lakes in Conflict Zones (50km)", "All National Lakes"))

ggplot(colombia_data, aes(x = species_richness, fill = Group)) +
  geom_density(alpha = 0.6, color = "white") +
  geom_rug(aes(color = Group), alpha = 0.3) +
  scale_fill_manual(values = c("All National Lakes" = "#FCD116",
                               "Lakes in Conflict Zones (50km)" = "#CE1126")) +
  scale_color_manual(values = c("All National Lakes" = "#FCD116",
                                "Lakes in Conflict Zones (50km)" = "#CE1126")) +
  labs(title = "Colombia",
       x = "Raw Species Richness (IUCN 2024)",
       y = "Density",
       fill = "Lake Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/colombia_species_richness_density.png", width = 10, height = 6)

# 6. Mexico -----
mexico_data <- lake_exposure_long %>%
  filter(year == 2024, country == "Mexico") %>%
  inner_join(lake_metadata, by = "Hylak_id") %>%
  mutate(Group = ifelse(is_exposed, "Lakes in Conflict Zones (50km)", "All National Lakes"))

ggplot(mexico_data, aes(x = species_richness, fill = Group)) +
  geom_density(alpha = 0.6, color = "white") +
  geom_rug(aes(color = Group), alpha = 0.3) +
  scale_fill_manual(values = c("All National Lakes" = "#006847",
                               "Lakes in Conflict Zones (50km)" = "#CE1126")) +
  scale_color_manual(values = c("All National Lakes" = "#006847",
                                "Lakes in Conflict Zones (50km)" = "#CE1126")) +
  labs(title = "Mexico",
       x = "Raw Species Richness (IUCN 2024)",
       y = "Density",
       fill = "Lake Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/mexico_species_richness_density.png", width = 10, height = 6)


# 7. Brazil ----
brazil_data <- lake_exposure_long %>%
  filter(year == 2024, country == "Brazil") %>%
  inner_join(lake_metadata, by = "Hylak_id") %>%
  mutate(Group = ifelse(is_exposed, "Lakes in Conflict Zones (50km)", "All National Lakes"))

# 2. Distribution Plot ----
ggplot(brazil_data, aes(x = species_richness, fill = Group)) +
  geom_density(alpha = 0.6, color = "white") +
  geom_rug(aes(color = Group), alpha = 0.3) +
  scale_fill_manual(values = c("All National Lakes" = "#009739",
                               "Lakes in Conflict Zones (50km)" = "#FEDD00")) +
  scale_color_manual(values = c("All National Lakes" = "#009739",
                                "Lakes in Conflict Zones (50km)" = "#FEDD00")) +
  labs(title = "Brazil",
       x = "Raw Species Richness (IUCN 2024)",
       y = "Density",
       fill = "Lake Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/brazil_species_richness_density.png", width = 10, height = 6)


