# Weighted Average Detection Efficiency (WADE) Positioning - Field Data Application
# Jake Brownscombe
# 04-09-2025

# this script uses positionR functions to apply WADE positioning to field data in a simple application.
# The main difference from the simulation example is there are helper functions to integrate commonly formatted (OTN)
# acoustic telemetry data into the WADE algorithm, and there is no ground truth (track) data to compare to.
# Users have options to integrate dynamic, site specific detection efficiency (range) models.


# 1. SETUP AND LIBRARIES =====================================================
library(positionR)
library(raster)
library(sf)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(patchwork)
library(lubridate)
library(tibble)


# Set seed for reproducibility
set.seed(123)


# 2. DATA LOADING ============================================================

# Load field data
data("stoney_fish_detections")  # Detection data from walleye in Stoney Lake
data("stoney_rx_deploy")        # Receiver deployment histories with coordinates
data("depth_raster")            # Bathymetry raster for study area
raster::crs(depth_raster) <- "EPSG:32617" #Set CRS for spatial data (WADE functions developed using UTM)

# 3. SPATIAL SETUP - CALCULATE DISTANCES AND DETECTION EFFICIENCY=============

# Calculate cost distances between all receiver stations and raster cells
# This creates the foundation for detection probability modeling
station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = stoney_rx_deploy,
  max_distance = 30000,  # Maximum distance to consider (m)
  station_col = "station_id"
)




# 4. DETECTION RANGE MODEL====================================================

# Detection range varies over space and time considerably; it's more ideal to use
# a more advanced detection range model from in situ sampling, which can be
# integrated into this workflow. Here we'll generate a simple model based on
# distance and depth, which is reasonably representative.

# Create logistic detection efficiency model based on distance and depth
logistic_DE <- create_logistic_curve_depth(
  min_depth = 1,
  max_depth = 35,
  d50_min_depth = 400,     # Distance at 50% detection efficiency in shallow water
  d95_min_depth = 800,     # Distance at 5% detection efficiency in shallow water
  d50_max_depth = 750,     # Distance at 50% detection efficiency in deep water
  d95_max_depth = 1500,    # Distance at 5% detection efficiency in deep water
  plot = TRUE,
  return_model = TRUE,
  return_object = TRUE
)

# Apply detection efficiency model to calculate probability for each receiver-cell combination
station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = station_distances %>%
    dplyr::rename(dist_m = cost_distance) %>%
    dplyr::mutate(depth_m = abs(raster_value)),
  type = "response"
)

# Visualize detection field for one receiver (example)
first_station <- stoney_rx_deploy$station_id[1]
ggplot(station_distances %>% filter(station_no == first_station),
       aes(x, y, fill = DE_pred)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Detection\nEfficiency") +
  geom_point(data = stoney_rx_deploy %>% filter(station_id == first_station),
             aes(x, y), col = "green", size = 1, inherit.aes = FALSE) +
  labs(title = paste("Detection Efficiency Field - Station", first_station),
       subtitle = "Green point indicates receiver location") +
  theme_minimal() +
  coord_sf()


# 5. SYSTEM-WIDE DETECTION EFFICIENCY=========================================


# Calculate cumulative detection probability across the entire receiver array
# This helps understand overall array coverage and can be used to generate absences within areas sampled
target_date <- as.POSIXct("2023-07-01", tz = "UTC")  # Select date (array changes over time)

system_DE <- calculate_detection_system(
  distance_frame = station_distances,
  receiver_frame = stoney_rx_deploy,
  model = logistic_DE$log_model,
  output_type = "cumulative",
  plots = TRUE,
  target_date = target_date
)


# 6. PREPARE DETECTION DATA FOR WADE POSITIONING==============================


# Select fish and time period for analysis
selected_fish <- "Walleye-1504321"  # Can be changed to analyze different fish
start_time <- as.POSIXct("2023-10-15 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2023-10-29 23:00:00", tz = "UTC")

# WADE requires a specific format of detections and non detections for stations deployed within time periods
wade_data <- prepare_detection_data_for_wade(
  fish_detections = stoney_fish_detections,
  station_deployments = stoney_rx_deploy,
  selected_fish_id = selected_fish,
  time_aggregation = "day",  # Daily aggregation for faster computation
  start_time = start_time,
  end_time = end_time
)

# Check for deployment warnings
if (length(wade_data$deployment_warnings) > 0) {
  cat("\n=== DEPLOYMENT WARNINGS ===\n")
  for (warning in wade_data$deployment_warnings) {
    cat(warning, "\n")
  }
  cat("\n")
}


# 7. WADE POSITIONING CALCULATIONS============================================

# Weighted Average Detection Efficiency (WADE) calculates a temporally aggregated positioning
#  estimate based on the weighted average (number of detections) of the detection efficiency
#fields, as well as non detection fields (receivers with no detections), and integrates these
# into a positioning estimate. The detection fields themselves may also be useful without
# considering non-detections.

# Apply WADE algorithm
positioning_results <- calculate_fish_positions(
  station_detections = wade_data$station_detections,
  station_distances_df = station_distances,
  station_info = wade_data$station_info,
  de_model = logistic_DE$log_model,
  detection_weight = 0.3,        # Weight for detected receivers
  non_detection_weight = 0.7,    # Weight for non-detected receivers. Overweighting non-dets as eg here
  max_non_detection_distance = 2000,  # Max distance to consider non-detections
  weighting_method = "normalize_stations",  # Method for combining probabilities
  normalization_method = "min_max", #normalization approach
  fish_id_col = "path_id",
  time_col = "datetime",
  time_aggregation = "day",
  station_col = "station_id",
  verbose = TRUE
)



# 8. VISUALIZATION OF POSITIONING RESULTS====================================

# Generate plots for available fish
head(positioning_results$position_probabilities %>% as.data.frame())
table(positioning_results$position_probabilities$time_period_posix)
cat("Available columns in position_probabilities:", paste(names(positioning_results$position_probabilities), collapse = ", "), "\n")
cat("Fish IDs in results:", unique(positioning_results$position_probabilities$fish_id), "\n")

# Select a time period with detections for visualization
available_times <- unique(positioning_results$position_probabilities$time_period_label)
selected_time <- available_times[12]  #selecting one eg time period
cat("Selected time period for plotting:", selected_time, "\n")

#plot positioning field (detection only)
plot_fish_positions(
        positioning_results = positioning_results,
        depth_raster_df = depth_raster,
        fish_select = selected_fish,  # Use original fish ID variable
        time_select = available_times[12],
        plot_type = "detection",
        prob_threshold = 0.01,  # Can change this to adjust which probability cells are displayed
        track_data = NULL  # No track data for field data - will plot probability fields only
      )

#plot positioning field (non detection)
plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = selected_fish,  # Use original fish ID variable
  time_select = available_times[12],
  plot_type = "non_detection",
  track_data = NULL
)

#plot positioning field (integrated)
plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = selected_fish,  # Use original fish ID variable
  time_select = available_times[12],
  plot_type = "integrated",
  prob_threshold = 0.7,
  detection_threshold = 0.01, #note to apply integrated, need to filter data down or else it will include other areas outside detection area**
  track_data = NULL)

plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = selected_fish,  # Use original fish ID variable
  time_select = available_times[13],
  plot_type = "detection",
  prob_threshold = 0.05,
  detection_threshold = 0.0, #note to apply integrated, need to filter data down or else it will include other areas outside detection area**
  track_data = NULL)

plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = selected_fish,  # Use original fish ID variable
  time_select = available_times[13],
  plot_type = "integrated",
  prob_threshold = 0.7,
  detection_threshold = 0.01, #note to apply integrated, need to filter data down or else it will include other areas outside detection area**
  track_data = NULL)

#more manual plotting style (all time periods here):

#prepare raster and detections objects
raster_df <- as.data.frame(depth_raster, xy = TRUE)
dets_sum <- wade_data$station_detections %>% group_by(path_id, datetime, station_id, station_x, station_y) %>% summarise(dets=sum(n_detections)) %>% rename(time_period_posix=datetime)
dets_sum$time_period_label <- as.character(dets_sum$time_period_posix)
dets_sum_det <- dets_sum %>% filter(dets>0)
dets_sum_nodet <- dets_sum %>% filter(dets==0)


ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  new_scale_fill() +
  geom_tile(data = positioning_results$position_probabilities %>%
              filter(weighted_mean_DE_normalized_scaled > 0.01 & integrated_prob>0.7),
            aes(x, y, fill = integrated_prob)) +
  scale_fill_viridis_c(option = "magma", name = "Probability", alpha = 0.8) +
  geom_point(data = dets_sum_det, aes(station_x, station_y, size = dets),
             col = "yellow", fill = NA, pch = 21) +
  geom_point(data = dets_sum_nodet, aes(station_x, station_y),
             col = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "WADE Estimates and Detections") +
  facet_wrap(~time_period_posix) +
  coord_sf()


ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  new_scale_fill() +
  geom_tile(data = positioning_results$position_probabilities %>%
              filter(weighted_mean_DE_normalized_scaled > 0.01 & integrated_prob>0.7 & time_period_label=="2023-10-27"),
            aes(x, y, fill = integrated_prob), alpha=1) +
  scale_fill_viridis_c(option = "magma", name = "Probability", alpha = 0.8) +
  geom_point(data = dets_sum_det %>% filter(time_period_label=="2023-10-27"),
             aes(station_x, station_y, size = dets), col = "yellow", fill = NA, pch = 21, stroke=2) +
  geom_point(data = dets_sum_nodet %>% filter(time_period_label=="2023-10-27"),
             aes(station_x, station_y), col = "red", size = 2) +
  theme_minimal() +
  labs(title = "WADE Estimates and Detections") +
  facet_wrap(~time_period_posix) +
  coord_sf(xlim=c(738000, 730000), ylim=c(4936000, 4942000))



# 9. SPACE USE ESTIMATION FROM WADE POSITIONS=================================


#pre filter data based on selected probability thresholds (flexible options here)
#also rescaling the positioning field (integrated_prob selected here) to 0-1 after filtering,
#which works better for generating a reasonable point distribution
filtered_results <- positioning_results
filtered_results$position_probabilities <- positioning_results$position_probabilities %>%
  filter(integrated_prob > 0.7,                        # integrated threshold
         weighted_mean_DE_normalized_scaled > 0.01) %>% # detection threshold
  mutate(
    # Rescale integrated_prob to use full [0,1] range after filtering
    integrated_prob_min = min(integrated_prob, na.rm = TRUE),
    integrated_prob_max = max(integrated_prob, na.rm = TRUE),
    integrated_prob_range = integrated_prob_max - integrated_prob_min,
    integrated_prob_rescaled = ifelse(integrated_prob_range > 0,
                                      (integrated_prob - integrated_prob_min) /
                                        integrated_prob_range,
                                      0.5)  # Fallback if no variation
  ) %>%
  select(-integrated_prob_min, -integrated_prob_max, -integrated_prob_range)

#generate points based on rescaled probability
position_points <- sample_points_from_probabilities(
  filtered_results,
  prob_column = "integrated_prob_rescaled",
  n_points = 1000,
  crs = 32617,
  by_group = TRUE,
  seed = 456
)

#plot
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = position_points, aes(color = probability),
          alpha = 0.5, size = 1) +
  scale_color_viridis_c(name = "Probability", option="magma") +
  geom_point(data=dets_sum_det, aes(station_x, station_y, size=dets), col="yellow", fill="NA", pch=21)+
  geom_point(data=dets_sum_nodet, aes(station_x, station_y), col="red", size=0.5)+
  theme_minimal()+labs(title="Presence Points")+facet_wrap(~time_period_posix)


#aggregate by location to double check we're generating points relative to the position prob:
# Extract coordinates and count points

position_points_count <- position_points %>%
  mutate(coords = st_coordinates(.),
         x = coords[,1],
         y = coords[,2]) %>%
  st_drop_geometry() %>%
  group_by(x, y, fish_id, time_period_posix, time_period_label) %>%
  summarise(point_count = n(), .groups = "drop")

# Create the plot with point counts as color
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  new_scale_fill() +
  geom_point(data = position_points_count, aes(x = x, y = y, color = point_count),
             alpha = 0.7, size = 1.2) +
  scale_color_viridis_c(name = "Point Count", option = "magma") +
  geom_point(data = dets_sum_det, aes(station_x, station_y, size = dets),
             col = "yellow", fill = NA, pch = 21) +
  geom_point(data = dets_sum_nodet, aes(station_x, station_y),
             col = "red", size = 0.5) +
  theme_minimal() +
  coord_sf()+
  labs(title = "Sampled Point Density") +
  facet_wrap(~time_period_posix)


#these points are directly useful for habitat selection (below).
#here we'll calculate 'home ranges' using constrained convex hull,
#to get a sense of the scale of space use:

# Extract coordinates from sf object
position_points_coords <- position_points %>%
  dplyr::mutate(coords = sf::st_coordinates(geometry),
                x = coords[,1], y = coords[,2]) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-coords)


# calculate space use with constrained convex hull
position_space_use_results <- calculate_space_use(
  track_data = position_points_coords,
  prob_column = "probability",
  prob_thresholds = 0.05, #can get a sense of this value from simulations
  by_fish = TRUE,
  by_time_period = TRUE,
  time_aggregation = "day",
  methods = c("constrained_convex_hull"),
  grid_resolution = 100,
  reference_raster = depth_raster
)


#plot
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  new_scale_fill() +
  geom_tile(data = position_space_use_results$spatial_data$constrained_convex_hull %>%
              rename(time_period_posix=time_period_label),
            aes(x, y), fill="green") +
  geom_point(data = dets_sum_det, aes(station_x, station_y, size = dets),
             col = "yellow", fill = NA, pch = 21) +
  geom_point(data = dets_sum_nodet, aes(station_x, station_y),
             col = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "WADE Home Range Estimates and Detections") +
  facet_wrap(~time_period_posix) +
  coord_sf()



# 10. HABITAT SELECTION FROM WADE POSITIONS===================================


#have presence points from above:
head(position_points)

#generate absences
#different options here, will sample from within the system DE uniformly with
#a cutoff detection probability threshold
absence_points_uniform <- sample_points_from_system_de(
  system_DE,
  position_points = position_points,
  n_points = 1000,
  uniform = TRUE,  #set uniform points
  min_prob_threshold = 0.1, #DE cutoff
  crs = 32617,
  seed = 123
)


#plot presence/absence points:
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = position_points  %>% filter(fish_id=="Walleye-1504321" & time_period_label==selected_time), aes(color = probability),
          alpha = 0.5, size = 1) +
  scale_color_viridis_c(name = "Probability", option="magma") +
  theme_minimal()+labs(title="Presence Points")+

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = absence_points_uniform %>% filter(fish_id=="Walleye-1504321" & time_period_label==selected_time),
          alpha = 1, size = 0.8, col="red") +
  scale_color_viridis_c(name = "Cumulative\nProbability") +
  theme_minimal()+labs(title="Absence Points (Uniform System DE Filtered)")



#depth selection
position_points$depth_m <- raster::extract(depth_raster, position_points)
absence_points_uniform$depth_m <- raster::extract(depth_raster, absence_points_uniform)


#Combine presence and absence data
presence_absence_points <- rbind(
  position_points %>%
    select(fish_id, time_period_posix, depth_m) %>%
    mutate(type = "presence"),
  absence_points_uniform %>%
    select(fish_id, time_period_posix, depth_m) %>%
    mutate(type = "absence")
)


plot_depth_selection(presence_absence_points, plot_type = "density") +
  ggtitle("Depth Selection") +
  labs(subtitle = "Blue = presences, Red = absences") +
  #by time
  plot_depth_comparison(presence_absence_points,
                        comparison_var = "time_period_posix",
                        plot_type = "boxplot") +
  ggtitle("Depth Selection by Time Period")


#at this stage you can fit a model of your choosing to assess habitat selection

#fin.

