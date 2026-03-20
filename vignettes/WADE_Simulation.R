# Weighted Average Detection Efficiency (WADE) Positioning on Simulation Data
# Jake Brownscombe
# 04-09-2025

# this script uses positionR functions to run Weighted Average Detection Efficiency (WADE)
# positioning on acoustic telemetry data. This example runs off of simulated data
# (arrays, detection models, movement tracks). It may be particularly useful for:
# - assessing array performance
# - determining optimal approaches with WADE for field data


# 1. SETUP AND LIBRARIES =====================================================
library(positionR)
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set seed for reproducibility
set.seed(123)

# 2. DATA LOADING ============================================================

# Load simulation environment data
data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617" # Set CRS for spatial data (WADE functions developed using UTM)

# 3. ARRAY DESIGN ============================================================

# Generate receiver array with regular spacing
points_regular <- generate_exact_regular_points(depth_raster, n_points = 80, seed = 123) # Regular spacing of set # of points (actual number may vary slightly)

# Visualize array design
plot_points_on_input(depth_raster, points_regular)

# Calculate cost distances between receivers and raster cells
# This creates the foundation for detection probability modeling
station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = points_regular,
  max_distance = 30000,  # Maximum distance to consider (30 km)
  station_col = "station_id"
)


# 4. DETECTION RANGE MODEL ===================================================

# Create depth-dependent detection efficiency model
# These parameters should be based on empirical range testing in your study system
logistic_DE <- create_logistic_curve_depth(
  min_depth = 1,           # Minimum depth (m)
  max_depth = 35,          # Maximum depth (m)
  d50_min_depth = 400,     # 50% detection range at shallow depth (m)
  d95_min_depth = 800,     # 5% detection range at shallow depth (m)
  d50_max_depth = 750,     # 50% detection range at deep depth (m)
  d95_max_depth = 1500,    # 5% detection range at deep depth (m)
  plot = TRUE,
  return_model = TRUE,
  return_object = TRUE
)

# Apply detection range model to predict detection efficiency for each receiver
# These predictions are essential for calculating animal positions and habitat selection
station_distances$DE_pred <- stats::predict(logistic_DE$log_model,
                                            newdata = station_distances %>%
                                            dplyr::rename(dist_m = cost_distance) %>%
                                            dplyr::mutate(depth_m = abs(raster_value)),
                                            type = "response")

# Visualize detection efficiency field for example station
ggplot(station_distances %>% dplyr::filter(station_no == 10),
       aes(x, y, fill = DE_pred)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  ggtitle("Detection Efficiency Field (Station 10)") +
  theme_minimal() +
  coord_sf()



# 5. SYSTEM DETECTION COVERAGE ===============================================

# Calculate cumulative system detection efficiency
# This will be used for generating absence points in habitat selection analysis
# include_barriers = TRUE incorporates barrier masking into coverage calculations
system_DE <- calculate_detection_system(
  distance_frame = station_distances,
  receiver_frame = points_regular,
  model = logistic_DE$log_model,
  output_type = "cumulative",  # Probability any receiver detects
  plots = TRUE,
  include_barriers = TRUE  # Prevent detections through land obstacles
)







# 6. FISH MOVEMENT SIMULATION ================================================

# Set start time to maintain temporal structure throughout analysis
start_time <- as.POSIXct("2025-07-15 12:00:00", tz = "UTC")

# Generate fish tracks with detections
# include_barriers = TRUE prevents detections through land obstacles
fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,  # Contains DE_pred values
  n_paths = 5,
  n_steps = 480,
  step_length_mean = 50,
  step_length_sd = 30,
  time_step = 180,
  seed = 123,
  start_time = start_time,
  include_barriers = TRUE  # Prevent detections through land
)



# Visualize fish tracks with detection events
plot_fish_tracks(
  fish_simulation,    # Simulation outputs
  depth_raster,       # Depth raster
  points_regular,     # Receiver array used
  show_detections = TRUE
)

# Analyze detection performance across the array
detection_performance <- analyze_detection_performance(fish_simulation,
                                                       create_plots = TRUE,
                                                       display_plots = FALSE)

# Display performance plots
detection_performance$plots$by_path +
  detection_performance$plots$by_station +
  detection_performance$plots$distribution +
  detection_performance$plots$time_series


# 7. WADE POSITION ESTIMATION ================================================

# Apply WADE algorithm to estimate fish positions
# Integrates both detection and non-detection evidence
# include_barriers = TRUE prevents position estimates through land obstacles,
# improving accuracy by accounting for physical geography
positioning_results <- calculate_fish_positions(
  station_detections = fish_simulation$station_detections,
  station_distances_df = station_distances,
  station_info = points_regular,
  de_model = logistic_DE$log_model,
  integration_method = "subtractive",  # "subtractive", "multiplicative", or "additive"
  max_non_detection_distance = 2000,   # Maximum range for non-detection inference (m)
  weighting_method = "normalize_stations",  # Normalize across stations
  # weighting_method = "information_theoretic",  # Alternative weighting approach
  normalization_method = "min_max",    # Min-max normalization
  fish_id_col = "path_id",
  time_col = "datetime",
  time_aggregation = "day",            # Daily position estimates
  station_col = "station_id",
  include_barriers = TRUE,             # Account for barriers in positioning
  verbose = TRUE
)



# Visualize position estimates for example fish and time period
fish_select <- 1
time_select <- "2025-07-15"

# Plot detection-based positioning
plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "detection",
  prob_threshold = 0.01,  # Display cells above 1% probability
  track_data = fish_simulation$tracks
) +

# Plot non-detection-based positioning
plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "non_detection",
  track_data = fish_simulation$tracks
) +

# Plot integrated positioning (combined evidence)
plot_fish_positions(
  positioning_results = positioning_results,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "integrated",
  prob_threshold = 0.4,
  detection_threshold = 0.01,  # Filter to detection area only
  track_data = fish_simulation$tracks
) +
  plot_layout(ncol = 2)



# 8. SPACE USE ANALYSIS =======================================================

#first generate points within the probability fields for each fish/time period
position_points <- sample_points_from_probabilities(
  positioning_results,
  prob_column = "weighted_mean_DE_normalized_scaled",
  n_points = 1000,
  min_prob_threshold = 0.0, #no prob cutoff here, will deal with this below
  crs = 32617,
  by_group = TRUE,
  seed = 456
)

# Extract coordinates from sf object
position_points_coords <- position_points %>%
  dplyr::mutate(coords = sf::st_coordinates(geometry),
    x = coords[,1], y = coords[,2]) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-coords)


# calculate space use
position_space_use_results <- calculate_space_use(
  track_data = position_points_coords,
  prob_column = "probability",
  prob_thresholds = 0.2, #can alter this to compare size of area to actual track, including multiple (see below)
  by_fish = TRUE,
  by_time_period = TRUE,
  time_aggregation = "day",
  methods = c("constrained_convex_hull"),
  grid_resolution = 100,
  reference_raster = depth_raster
)



#actual track area
track_space_use_results <- calculate_space_use(
  track_data = fish_simulation$tracks,
  by_fish = TRUE,
  by_time_period = TRUE,
  time_aggregation = "day",  # Options: "hour", "day", "month", "none"
  methods = c("grid_cell_count", "constrained_convex_hull"),
  grid_resolution = 100,
  reference_raster = depth_raster
)


#plot
fish_select <- 1
time_select <- "2025-07-15"
plot_space_use(position_space_use_results, #not working - the same as track space use**??
               plot_type = "map",
               fish_select = fish_select,
               time_select = time_select,
               method_select = "constrained_convex_hull",
               track_data = fish_simulation$tracks,
               background_raster = depth_raster) +
  ggtitle("Constrained Convex Hull Space Use (Positioning)")+

plot_space_use(track_space_use_results,
               plot_type = "map",
               track_data = fish_simulation$tracks,
               fish_select = fish_select,
               time_select = time_select,
               method_select = "constrained_convex_hull",
               background_raster = depth_raster) +
  ggtitle("Constrained Convex Hull Space Use (Track)")



#compare space use and probability thresholds
multi_results <- calculate_space_use(
  track_data = position_points_coords,
  prob_column = "probability",
  by_fish = TRUE,
  by_time_period = TRUE,
  prob_thresholds = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
  methods = c("constrained_convex_hull"),
  reference_raster = depth_raster
)

# Compare against track-based space use
comparison <- compare_space_use_thresholds( #ISSUE here, track overlaps all come out as 0.
  multi_threshold_results = multi_results,
  reference_tracks = fish_simulation$tracks,
  reference_space_use = track_space_use_results
)


# View results
print(comparison$plots$area_ratio_boxplot) +
print(comparison$plots$track_overlap_boxplot)
#0.4 appears suitable for space use area (home range), covers ~80% of the actual track.


#calculate space use with this threshold
space_use_0.4 <- calculate_space_use(
  track_data = position_points_coords,
  prob_column = "probability",
  by_fish = TRUE,
  by_time_period = TRUE,
  prob_thresholds =  0.4,
  methods = c("constrained_convex_hull"),
  reference_raster = depth_raster
)

ggplot(space_use_0.4$space_use_estimates, aes(as.factor(fish_id), constrained_convex_hull_area_hectares))+
  geom_boxplot()+theme_minimal()+labs(title="Positioning Estimates of Space Use", x="Fish ID", y="Space Use (hectares)")




# 9. HABITAT SELECTION ANALYSIS ==============================================

# Generate presence points for habitat selection analysis
# Note: For integrated_prob values, pre-filter based on detection threshold
# to ensure areas are within the positioning field. Here we use the detection field
# (weighted_mean_DE_normalized_scaled) with optimized threshold:
position_points <- sample_points_from_probabilities(
  positioning_results,
  prob_column = "weighted_mean_DE_normalized_scaled",
  n_points = 1000,
  min_prob_threshold = 0.4, #based on overlap results above
  crs = 32617,
  by_group = TRUE,
  seed = 456
)


#plot one eg
raster_df <- as.data.frame(depth_raster, xy = TRUE)
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = position_points %>% filter(fish_id==1 & time_period_label=="2025-07-15"), aes(color = probability),
          alpha = 0.5, size = 1) +
  scale_color_viridis_c(name = "Probability", option="magma") +
  theme_minimal()+labs(title="Presence Points")


#generate absences
#here we'll use the positioning system area as a measure of what is actually being sampled
absence_points <- sample_points_from_system_de(
  system_DE,
  position_points = position_points,
  n_points = 1000,
  min_prob_threshold = 0.05,
  crs = 32617,
  seed = 123
)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = absence_points %>% filter(fish_id==1 & time_period_label=="2025-07-15"), aes(color = probability),
          alpha = 1, size = 0.8) +
  scale_color_viridis_c(name = "Cumulative\nProbability") +
  theme_minimal()+labs(title="Absence Points (System DE relative)")


#or can sample uniformly across the whole raster:
absence_points_uniform <- sample_points_from_system_de(
  system_DE,
  position_points = position_points,
  n_points = 1000,
  uniform = TRUE,  #set uniform points
  min_prob_threshold = 0.0, #zero threshold, whole system
  crs = 32617,
  seed = 123
)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = absence_points_uniform %>% filter(fish_id==1 & time_period_label=="2025-07-15"),
          alpha = 1, size = 0.8, col="red") +
  scale_color_viridis_c(name = "Cumulative\nProbability") +
  theme_minimal()+labs(title="Absence Points (Uniform Across System)")




#depth selection analysis
position_points$depth_m <- raster::extract(depth_raster, position_points)
absence_points_uniform$depth_m <- raster::extract(depth_raster, absence_points_uniform)


#Combine presence and absence data
#will use uniform (whole system) absences here, mainly for consistency when
#comparing habitat selection with the actual track below
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
  labs(subtitle = "Blue = presences, Red = absences")+

  # Compare depth selection by fish
  plot_depth_comparison(presence_absence_points,
                        comparison_var = "fish_id",
                        plot_type = "boxplot") +
  ggtitle("Depth Selection by Individual Fish")








# 10. TRACK-BASED HABITAT SELECTION (GROUND TRUTH) ==========================

# Generate presence points from actual track data (ground truth)
track_presence <- sample_points_from_track_grid(
  track_data = fish_simulation$tracks,
  reference_raster = depth_raster,
  n_points = 1000,
  crs = 32617
)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = track_presence, aes(color = count),
          alpha = 0.5, size = 1) +
  scale_color_viridis_c(name = "Count", option = "magma") +
  coord_sf() +
  theme_minimal() +
  labs(title = "Track Presence Points")


#track absences
track_absence<- sample_points_from_system_de(
  system_DE,
  position_points = track_presence,
  n_points = 1000,
  uniform = TRUE,  #set uniform points across raster
  min_prob_threshold = 0.0, #zero threshold, whole system
  crs = 32617,
  seed = 123
)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_sf(data = track_absence %>% filter(fish_id==1 & time_period_label=="2025-07-15"),
          alpha = 1, size = 0.8, col="red") +
  scale_color_viridis_c(name = "Cumulative\nProbability") +
  theme_minimal()+labs(title="Track Absence Points (Uniform Across System)")


#track depth selection
track_presence$depth_m <- raster::extract(depth_raster, track_presence)
track_absence$depth_m <- raster::extract(depth_raster, track_absence)

# Combine presence and absence data
track_presence_absence_points <- rbind(
  track_presence %>%
    select(fish_id, time_period_posix, depth_m) %>%
    mutate(type = "presence"),
  track_absence %>%
    select(fish_id, time_period_posix, depth_m) %>%
    mutate(type = "absence")
)

plot_depth_selection(track_presence_absence_points, plot_type = "density") +
  ggtitle("Depth Selection") +
  labs(subtitle = "Blue = presences, Red = absences")+

  # Compare depth selection by fish
  plot_depth_comparison(track_presence_absence_points,
                        comparison_var = "fish_id",
                        plot_type = "boxplot") +
  ggtitle("Depth Selection by Individual Fish")




# 11. COMPARATIVE HABITAT SELECTION ANALYSIS ================================

# Compare habitat selection patterns between positioning estimates and ground truth
# Uses Random Forest models with partial dependence plots

comparative_results <- analyze_comparative_habitat_selection(
  positioning_data = presence_absence_points,
  space_use_data = track_presence_absence_points,
  analysis_type = "both",
  formula = presence ~ depth_m + fish_id + time_period_posix,
  time_variable = "time_period_posix",
  sample_size = 10000,
  ntree = 500,
  create_plots = TRUE,
  create_comparison = TRUE
)

# View comprehensive summary
print_comparative_summary(comparative_results)


# Plot all results (individual + comparison)
plot_comparative_results(comparative_results, plot_type = "comparison")


#With optimal thresholds, can estimate both home range and habitat selection effectively (88%) with this
#underlying movement data, receiver array, and DE model

#fin.
