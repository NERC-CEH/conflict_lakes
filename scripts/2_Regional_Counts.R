#Script 2:
#  • Join the exposure results from script 1 with country and ipbes regional polygons to calc % exposed lakes per region per year
#  • If a county is at war, how many additonal lakes (not in direct exposure zones) are there?


# 1. Libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(FAOSTAT)
library(httr)
library(graticule)

# 2. Load and Tag Data ----
load("data/globo_topo_poly.rda")
conflict <- read_csv("data/GEDEvent_v25_1.csv")

# read in all hydrolake points (huge)
hydro_gdb <- "data/HydroLAKES_points_v10_shp"
st_layers(hydro_gdb)
lake_centroids <- st_read(hydro_gdb, layer = "HydroLAKES_points_v10")

# create world tags to look at global regions
world_tags <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(name_long, region_un, continent) %>%
  st_transform(54009)

# calcutlate lake centroids
lake_centroids_tagged <- lake_centroids %>%
  st_make_valid() %>%
  #st_centroid() %>%
  st_transform(54009) %>%
  st_join(world_tags) %>%
  rename(country = name_long, region = region_un)

# 3. Generate Annual Exposure   ----
target_buffer <- 50000 # 50km
anchor_years <- 2000:2024 # contrinutous 2000-24

conflict_sf <- conflict %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(54009)

message("Generating annual exposure key... this may take a moment.")
lake_exposure_long <- map_df(anchor_years, function(yr) {
  # filter conflict buffers for each year
  conf_yr_buffer <- conflict_sf %>%
    filter(year == yr) %>%
    st_buffer(target_buffer) %>%
    st_union()

  #  which lake centroids are inside the buffer
  is_exposed <- st_intersects(lake_centroids_tagged, conf_yr_buffer, sparse = FALSE)

  lake_centroids_tagged %>%
    st_drop_geometry() %>%
    select(Hylak_id, country, region) %>%
    mutate(year = yr, is_exposed = as.vector(is_exposed))
})

# 4. Regional Exposure to Conflict ----
#  % of lakes exposed per region and per year.
regional_exposure <- lake_exposure_long %>%
  group_by(region, year) %>%
  summarise(
    total_lakes = n(),
    exposed_lakes = sum(is_exposed),
    pct_exposed = (exposed_lakes / total_lakes) * 100,
    .groups = "drop")

# 5. Non-exposed lakes in countries 'at war' ----
#  >= 25 d deaths (UCDP standard def for war)
war_countries <- conflict %>%
  group_by(year, country) %>%
  summarise(total_deaths = sum(best), .groups = "drop") %>%
  filter(total_deaths >= 25)

# seperate direct vs indirect exposure categories
additional_lakes_analysis <- lake_exposure_long %>%
  inner_join(war_countries, by = c("year", "country")) %>%
  mutate(
    status = case_when(
      is_exposed == TRUE ~ "Direct Exposure (within 50km)",
      is_exposed == FALSE ~ "Indirect Exposure (Country at War)"))

# 6. Plots  ----
war_lake_counts <- additional_lakes_analysis %>%
  group_by(year, status) %>%
  summarise(lake_count = n(), .groups = "drop")

countries_war <- ggplot(war_lake_counts, aes(x = year, y = lake_count, fill = status)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("Direct Exposure (within 50km)" = "#d73027",
                               "Indirect Exposure (Country at War)" = "#fc8d59")) +
  labs(title = "Lakes in Conflict-Affected Countries (2000-2024)",
       y = "Number of Lakes", x = "Year", fill = "Exposure Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ggsave
ggsave("output/lakes_in_war_countries.png", countries_war, width = 10, height = 6)

# same plot but ONLY country = ukraine
ukraine_lakes <- additional_lakes_analysis %>%
  filter(country == "Ukraine") %>%
  group_by(year, status) %>%
  summarise(lake_count = n(), .groups = "drop")

(ukraine_plot <- ggplot(ukraine_lakes, aes(x = year, y = lake_count, fill = status)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("Direct Exposure (within 50km)" = "#d73027",
                               "Indirect Exposure (Country at War)" = "#fc8d59")) +
  labs(title = "Lakes in Ukraine (2000-2024)",
       #subtitle = "Distinguishing between direct physical proximity and broader national instability",
       y = "Number of Lakes", x = "Year", fill = "Exposure Type") +
  theme_minimal() +
  theme(legend.position = "bottom"))

ggsave("output/ukraine_lakes.png", ukraine_plot, width = 10, height = 6)


# try regions ----
world_tags <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(name_long, region_un, subregion, continent) %>% # Added subregion
  st_transform(54009)

lake_centroids_tagged <- globo_topo_poly %>%
  st_make_valid() %>%
  st_centroid() %>%
  st_transform(54009) %>%
  st_join(world_tags)

na_lakes <- is.na(lake_centroids_tagged$region_un)

if(any(na_lakes)) {
  nearest_indices <- st_nearest_feature(lake_centroids_tagged[na_lakes, ], world_tags)
  lake_centroids_tagged[na_lakes, c("name_long", "region_un", "subregion", "continent")] <-
    st_drop_geometry(world_tags[nearest_indices, ])
}

lake_centroids_tagged <- lake_centroids_tagged %>%
  rename(country = name_long, region = region_un)
# plot of regional % lakes exposed to conflict
(regional_exposure_plot <- regional_exposure %>%
    filter(!is.na(region)) %>% # Safety filter
    ggplot(aes(x = year, y = pct_exposed, color = region)) + # Using 'region' as created above
    geom_line(size = 1) +
    labs(title = "Regional Lake Exposure to Conflict (2000-2024)",
         x = "Year", y = "% of Lakes Exposed in 50km buffer", color = "Region") +
    theme_minimal() +
    theme(legend.position = "bottom"))

# ggsave
ggsave("output/regional_exposure_plot.png", regional_exposure_plot, width =
       10, height = 6)




















