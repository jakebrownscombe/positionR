# Goal-directed movement simulations (biased correlated random walk)
# Jake Brownscombe
# 2026-03-12

# This vignette demonstrates goal-directed movement using a biased correlated random walk (BCRW).
# Fish navigate from a start location toward a goal with tunable bias strength (goal_bias).
# A bias of 0 produces a pure correlated random walk, while 1 drives the fish in a straight line.
# The fish stops moving when it arrives within goal_tolerance of the destination.


library(positionR)
library(raster)
library(dplyr)
library(ggplot2)


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


# --- Pick start and goal coordinates ---
# Extract valid water cells and choose two well-separated points
water_cells <- raster::as.data.frame(depth_raster, xy = TRUE)
water_cells <- water_cells[!is.na(water_cells[, 3]), ]

# Pick points near opposite ends of the lake (min and max x-coordinates)
start_idx <- which.min(water_cells$x)
goal_idx <- which.max(water_cells$x)

start_loc <- matrix(c(water_cells$x[start_idx], water_cells$y[start_idx]), ncol = 2)
goal_loc <- matrix(c(water_cells$x[goal_idx], water_cells$y[goal_idx]), ncol = 2)

cat("Start:", start_loc[1, 1], ",", start_loc[1, 2], "\n")
cat("Goal: ", goal_loc[1, 1], ",", goal_loc[1, 2], "\n")
cat("Euclidean distance:", round(sqrt(sum((start_loc - goal_loc)^2)) / 1000, 1), "km\n")



# ============================================================================
# Example 1: Default goal-directed walk (goal_bias = 0.5)
# ============================================================================

sim1 <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,
  n_paths = 1,
  n_steps = 5000,
  time_step = 60,
  start_locations = start_loc,
  goal_locations = goal_loc,
  goal_bias = 0.5,
  seed = 42
)

# Check endpoint relative to goal
tail(sim1$tracks, 3)

# Plot the track with start/goal markers
plot_fish_tracks(sim1,
                 depth_raster,
                 stations_regular,
                 show_detections = TRUE) +
  geom_point(aes(x = start_loc[1, 1], y = start_loc[1, 2]), shape = 21, fill = "green", size = 3) +
  geom_point(aes(x = goal_loc[1, 1], y = goal_loc[1, 2]), shape = 21, fill = "red", size = 3)



# ============================================================================
# Example 2: Bias strength comparison (0, 0.25, 0.5, 0.75, 1.0)
# ============================================================================

bias_values <- c(0, 0.25, 0.5, 0.75, 1.0)
all_tracks <- list()

for (i in seq_along(bias_values)) {
  sim_i <- simulate_fish_tracks(
    raster = depth_raster,
    station_distances = station_distances,
    n_paths = 1,
    n_steps = 5000,
    time_step = 60,
    start_locations = start_loc,
    goal_locations = goal_loc,
    goal_bias = bias_values[i],
    seed = 42
  )
  tracks_i <- sim_i$tracks
  tracks_i$bias <- paste0("bias = ", bias_values[i])
  all_tracks[[i]] <- tracks_i
}

bias_tracks <- do.call(rbind, all_tracks)
bias_tracks$bias <- factor(bias_tracks$bias, levels = paste0("bias = ", bias_values))

# Custom ggplot: simple background + coloured paths by bias level
library(ggplot2)
raster_df <- raster::as.data.frame(depth_raster, xy = TRUE)
names(raster_df)[3] <- "depth"
raster_df <- raster_df[!is.na(raster_df$depth), ]

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y), fill = "grey85") +
  geom_path(data = bias_tracks, aes(x = x, y = y, colour = bias), linewidth = 0.6, alpha = 0.8) +
  geom_point(aes(x = start_loc[1, 1], y = start_loc[1, 2]), shape = 21, fill = "green", size = 3) +
  geom_point(aes(x = goal_loc[1, 1], y = goal_loc[1, 2]), shape = 21, fill = "red", size = 3) +
  scale_colour_viridis_d(name = "Goal bias") +
  labs(title = "Effect of goal_bias on path directedness",
       x = "Easting (m)", y = "Northing (m)") +
  coord_equal() +
  theme_minimal()


# ============================================================================
# Example 3: Multi-path goal-directed simulation
# ============================================================================
# 5 paths, all sharing the same start and goal — shows natural path variation

n <- 5
start_locs <- matrix(rep(start_loc, each = n), ncol = 2)
goal_locs <- matrix(rep(goal_loc, each = n), ncol = 2)

sim3 <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,
  n_paths = n,
  n_steps = 5000,
  time_step = 60,
  start_locations = start_locs,
  goal_locations = goal_locs,
  goal_bias = 0.5,
  seed = 42
)

ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y), fill = "grey85") +
  geom_path(data = sim3$tracks, aes(x = x, y = y, colour = factor(path_id)), linewidth = 0.6, alpha = 0.8) +
  geom_point(aes(x = start_loc[1, 1], y = start_loc[1, 2]), shape = 21, fill = "green", size = 3) +
  geom_point(aes(x = goal_loc[1, 1], y = goal_loc[1, 2]), shape = 21, fill = "red", size = 3) +
  scale_colour_viridis_d(name = "Path") +
  labs(title = "Multi-path goal-directed simulation (goal_bias = 0.5)",
       x = "Easting (m)", y = "Northing (m)") +
  coord_equal() +
  theme_minimal()
#fin.

