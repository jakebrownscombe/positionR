# Array Design and Movement Detection Simulation
# Jake Brownscombe
# 04-09-2025

#this script uses positionR functions to:
# - generate receiver arrays
# - generate basic detection efficiency models
# - generate estimates of detection efficiency across arrays
# - run movement simulations and generate detection patterns
# - conduct space use analyses
# - generate presence absence data for habitat selection analyses



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

# Load study area environmental data
data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617" # Set CRS for spatial data (UTM Zone 17N required for distance calculations)


# 3. RECEIVER ARRAY DESIGN ===================================================

# Regular spacing - specified number of points in systematic pattern
stations_regular <- generate_regular_points(depth_raster, n_points = 50, seed = 123)

# Fixed spacing - points separated by specified distance
stations_spaced <- generate_spaced_points(depth_raster, spacing = 1000, seed = 123)

# Random locations - randomly distributed points
stations_random <- generate_random_points(depth_raster, n_points = 100, seed = 123)


# Visualize different array configurations
plot_points_on_input(depth_raster, stations_regular) +
  ggtitle("Regular Array Configuration")+

plot_points_on_input(depth_raster, stations_spaced) +
  ggtitle("Spaced Array Configuration (1000m)")

plot_points_on_input(depth_raster, stations_random) +
  ggtitle("Random Array Configuration")


# 4. DISTANCE FIELD CALCULATIONS ============================================

# Calculate cost-distance and straight-line distance from each receiver to all raster cells
# This creates the foundation for detection efficiency modeling across the study area
station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = stations_regular,  # Using regular array for analysis
  station_col = "station_id",
  max_distance = 30000  # 30 km maximum distance
)

# Visualize distance field for one receiver
ggplot(station_distances %>% filter(station_no == 1),
       aes(x, y, fill = cost_distance)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Distance (m)") +
  geom_point(data = stations_regular %>% filter(station_id == 1),
             aes(x, y), col = "green", size = 4, inherit.aes = FALSE) +
  labs(title = "Cost Distance Field from Station 1",
       subtitle = "Green point = station location") +
  theme_minimal() +
  coord_sf()

# Tortuosity visualization
ggplot(station_distances %>% filter(station_no == 1),
       aes(x, y, fill = tortuosity)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Tortuosity", limits = c(0.8, 1.6)) +
  geom_point(data = stations_regular %>% filter(station_id == 1),
             aes(x, y), col = "green", size = 4, inherit.aes = FALSE) +
  labs(title = "Tortuosity Field from Station 1",
       subtitle = "Green point = station location") +
  theme_minimal() +
  coord_sf()

# Barrier visualization
# The crosses_barrier column identifies where line-of-sight crosses land
# This prevents unrealistic detections through islands/peninsulas
ggplot(station_distances %>% filter(station_no == 1),
       aes(x, y, fill = crosses_barrier)) +
  geom_raster() +
  scale_fill_manual(values = c("#4A90A4", "#D4816F"),
                   labels = c("Clear path", "Crosses barrier")) +
  geom_point(data = stations_regular %>% filter(station_id == 1),
             aes(x, y), col = "green", size = 4, inherit.aes = FALSE) +
  labs(title = "Barrier Field from Station 1",
       subtitle = "Locations where line-of-sight crosses land") +
  theme_minimal() +
  coord_sf()


# 5. DETECTION RANGE MODEL ===================================================

# Create depth-dependent detection efficiency model
# In practice, use range testing data from your study system
# These parameters represent typical inland lake conditions
logistic_DE <- create_logistic_curve_depth(
  min_depth = 1,           # Minimum depth (m)
  max_depth = 35,          # Maximum depth (m)
  d50_min_depth = 400,     # 50% detection range at min depth (m)
  d95_min_depth = 800,     # 5% detection range at min depth (m)
  d50_max_depth = 750,     # 50% detection range at max depth (m)
  d95_max_depth = 1500,    # 5% detection range at max depth (m)
  plot = TRUE,
  return_model = TRUE,
  return_object = TRUE
)


# 6. DETECTION EFFICIENCY PREDICTION =========================================

# Apply the detection range model to predict detection efficiency
# for each receiver across the study area
station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = station_distances %>%
    rename(dist_m = cost_distance) %>%
    mutate(depth_m = abs(raster_value)),
  type = "response"
)

# Visualize detection efficiency field for one station
ggplot(station_distances %>% filter(station_no == 10),
       aes(x, y, fill = DE_pred)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Detection\nEfficiency") +
  labs(title = "Detection Efficiency Field - Station 10",
       subtitle = "Probability of detecting a transmission") +
  theme_minimal() +
  coord_sf()



# 7. SYSTEM DETECTION COVERAGE ===============================================

# Analyze detection performance at the system level to evaluate array design effectiveness
# include_barriers = TRUE incorporates barrier masking into system coverage calculations
system_DE <- calculate_detection_system(
  distance_frame = station_distances,
  receiver_frame = stations_regular,
  model = logistic_DE$log_model,
  output_type = "cumulative",  # Cumulative probability (any receiver detects)
  plots = TRUE,
  include_barriers = TRUE  # Prevent detections through land obstacles
)



# 8. FISH MOVEMENT SIMULATION ================================================

# Set start time for temporal consistency
start_time <- as.POSIXct("2025-07-15 12:00:00", tz = "UTC")

# Generate realistic fish movements using correlated random walks
# and simulate detections based on the array and detection model
# include_barriers = TRUE prevents detections through land obstacles,
# creating more realistic detection patterns
fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,  # Contains DE_pred values
  n_paths = 5,                          # Number of fish
  n_steps = 1000,                         # Steps per track
  step_length_mean = 50,                 # Average step length (m)
  step_length_sd = 30,                   # Step length variability
  time_step = 60,                        # Time between steps (seconds)
  seed = 1987,
  start_time = start_time,
  include_barriers = TRUE                # Prevent detections through land

  # Optional species-specific parameters:
  # mode = "species_empirical",    # or "species_theory" for 3-state behavioral model
  # species = "Walleye",           # Available: "Walleye", "Smallmouth Bass", "Muskellunge", "Generic"
)


# Visualize simulated tracks
raster_df <- as.data.frame(depth_raster, xy = TRUE)
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_path(data = fish_simulation$tracks,
            aes(x = x, y = y, group = path_id, color = factor(path_id)),
            linewidth = 0.8) +
  scale_color_discrete(name = "Fish ID") +
  labs(title = "Simulated Fish Movements",
       subtitle = "Correlated random walks in study area") +
  theme_minimal() +
  coord_sf()

# 9. MOVEMENT AND DETECTION ANALYSIS =========================================

# Plot tracks with detection information
plot_fish_tracks(
  fish_simulation,    # Simulation outputs
  depth_raster,       # Depth raster
  stations_regular,   # Receiver array
  show_detections = TRUE
) +labs(subtitle= "yellow circle = detected | red x = not detected")


# Analyze detection performance across the array
detection_performance <- analyze_detection_performance(fish_simulation)
# Display performance plots
detection_performance$plots$by_path +
  detection_performance$plots$by_station +
  detection_performance$plots$distribution +
  detection_performance$plots$time_series


# 10. SPACE USE ANALYSIS ======================================================

# Calculate space use estimates from the simulated tracks using multiple methods
space_use_results <- calculate_space_use(
  track_data = fish_simulation$tracks,
  by_fish = TRUE,
  by_time_period = TRUE,
  time_aggregation = "day",  # Options: "hour", "day", "month", "none"
  methods = c("convex_hull", "bounding_box", "grid_cell_count",
              "ellipse_95", "constrained_convex_hull", "mcp"),
  grid_resolution = 100,
  reference_raster = depth_raster
)


# Compare methods graphically
plot_space_use(space_use_results, plot_type = "comparison") +
  ggtitle("Space Use Method Comparison")


# Plot specific methods for selected fish and time
fish_select <- 1
time_select <- "2025-07-15"

# Convex hull method
plot_space_use(space_use_results,
               plot_type = "map",
               track_data = fish_simulation$tracks,
               fish_select = fish_select,
               time_select = time_select,
               method_select = "convex_hull",
               background_raster = depth_raster) +
  ggtitle("Convex Hull Space Use") +

# Constrained convex hull method (raster-based)
plot_space_use(space_use_results,
               plot_type = "map",
               track_data = fish_simulation$tracks,
               fish_select = fish_select,
               time_select = time_select,
               method_select = "constrained_convex_hull",
               background_raster = depth_raster) +
  ggtitle("Constrained Convex Hull Space Use") +

# Grid cell count method
plot_space_use(space_use_results,
               plot_type = "map",
               track_data = fish_simulation$tracks,
               fish_select = fish_select,
               time_select = time_select,
               method_select = "grid_cell_count",
               background_raster = depth_raster) +
  ggtitle("Grid Cell Count Space Use")



# 11. HABITAT SELECTION DATA GENERATION ======================================

# Generate presence points from track data for habitat selection analysis
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



# Generate absence points
track_absence<- sample_points_from_system_de(
  system_DE,
  position_points = track_presence,
  n_points = 1000,
  uniform = TRUE,  #set uniform points (not relative to system DE)
  min_prob_threshold = 0.1, #threshold cutoff - only points where >0.1 system DE
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





# 12. HABITAT SELECTION ANALYSIS =============================================

# Extract environmental values and analyze habitat selection patterns
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

#fin.
