# Advanced species movement simulations
# Jake Brownscombe
# 04-09-2025

# This is a beta version of this, with plans for expansion. The intention is to provide
# some basic mechanics and examples of how more advanced movement models can be integrated
# into movement simulations considering species and environmentally specific characteristics.
# there are some diverse options in how species can be parameterized to move into different
# behavioural states based on time of day and temperature, including spawning. The more advanced
# functionality relies on species_behavioral_params and daily_temperature files in the data folder.


library(positionR)


#load and process necessary data
data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617" # Set CRS for spatial data (UTM Zone 17N required for distance calculations)
stations_regular <- generate_exact_regular_points(depth_raster, n_points = 50, seed = 123)

#calc distances
station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = stations_regular,  # Using regular array for analysis
  station_col = "station_id",
  max_distance = 30000  # 30 km maximum distance
)

#DE model
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

#predict
station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = station_distances %>%
    rename(dist_m = cost_distance) %>%
    mutate(depth_m = abs(raster_value)),
  type = "response"
)



#run simulation
start_time <- as.POSIXct("2025-07-15 12:00:00", tz = "UTC")
walleye_sim <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,
  n_paths = 1,
  n_steps = 1000,
  time_step = 60,
  seed = 123,
  species = "Muskellunge",          # Available: "Walleye", "Smallmouth Bass", "Muskellunge"
  fish_size_cm = 45,                # Fish length in cm (defaults to typical_size_cm if not specified)
  behavioral_states = TRUE,         #three state movement model
  start_time = start_time,
  temperature_data = temporal_info
)

plot_fish_tracks(walleye_sim, #sim outputs
                 depth_raster, #depth raster
                 stations_regular, #receiver array used
                 show_detections = TRUE)


# Plot behavioral states with temperature
plot_behavioral_states(walleye_sim)


# Basic analysis with 5 temperature bins
temp_analysis <- analyze_behavioral_temperature(walleye_sim, by_fish_id = TRUE)

# View the plots
temp_analysis$plots$time_in_state
temp_analysis$plots$transition_probs



#spawning behaviour integrated
# Enable spawning behavior (requires species specification)
data("temporal_info")
start_time <- as.POSIXct("2025-02-01 12:00:00", tz = "UTC")

walleye_spawn_sim <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,
  n_paths = 1,
  n_steps = 30000,
  time_step = 120,
  seed = 123,
  species = "Walleye",
  spawning_behavior = TRUE,  # Enable spawning behavior
  temperature_data = temporal_info,  # Required for spawning
  start_time = start_time
)

# Check columns
names(walleye_spawn_sim$tracks)

#plots
plot_fish_tracks(walleye_spawn_sim, #sim outputs
                 depth_raster, #depth raster
                 show_detections = FALSE)

plot_fish_tracks(walleye_spawn_sim, #sim outputs
                 depth_raster, #depth raster
                 show_detections = FALSE,
                 color_by = "behavioral_state",
                 point_size = 0.3,
                 sample_rate = 1)

plot_fish_tracks(walleye_spawn_sim, #sim outputs
                 depth_raster, #depth raster
                 show_detections = FALSE,
                 color_by = "spawning_phase",
                 point_size = 0.3,
                 sample_rate = 1)

#fin.
