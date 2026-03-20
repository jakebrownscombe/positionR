# Particle Filter Positioning - Simulation Test
# Jake Brownscombe
# 2025

# This script tests the particle filter positioning algorithm using simulated fish tracks,
# comparing position estimates against ground truth and WADE estimates.

# 1. SETUP =====================================================================
library(positionR)
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)


# 2. ENVIRONMENT SETUP =========================================================

data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617"

# Generate receiver array
points_regular <- generate_exact_regular_points(depth_raster, n_points = 40, seed = 123)

# Calculate station distances
station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = points_regular,
  max_distance = 30000,
  station_col = "station_id"
)

# Create DE model
logistic_DE <- create_logistic_curve_depth(
  min_depth = 1, max_depth = 35,
  d50_min_depth = 400, d95_min_depth = 800,
  d50_max_depth = 750, d95_max_depth = 1500,
  plot = FALSE, return_model = TRUE, return_object = TRUE
)

# Apply DE model
station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = station_distances %>%
    dplyr::rename(dist_m = cost_distance) %>%
    dplyr::mutate(depth_m = abs(raster_value)),
  type = "response"
)


# 3. SIMULATE FISH TRACKS =====================================================

start_time <- as.POSIXct("2025-07-15 12:00:00", tz = "UTC")

fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,
  n_paths = 3,
  n_steps = 200,
  step_length_mean = 50,
  step_length_sd = 30,
  time_step = 180,  # 3-minute intervals
  seed = 42,
  start_time = start_time,
  include_barriers = TRUE
)


# 4. PARTICLE FILTER POSITIONING ===============================================

cat("\n=== Running Particle Filter ===\n")
pf_time <- system.time({
  pf_results <- particle_filter_positioning(
    detection_data = fish_simulation$station_detections,
    station_info = points_regular,
    de_model = logistic_DE$log_model,
    raster = depth_raster,
    n_particles = 1000,
    step_length_mean = 50,
    step_length_sd = 30,
    turning_angle_sd = 45,
    time_step = 180,
    max_distance = 30000,
    fish_id_col = "path_id",
    time_col = "datetime",
    station_col = "station_id",
    return_particles = FALSE,
    verbose = TRUE
  )
})
cat("Particle filter time:", pf_time["elapsed"], "s\n")
cat("Position estimates:", nrow(pf_results$positions), "\n")


# 5. COMPARE WITH GROUND TRUTH ================================================

# Get true positions from simulation tracks
true_tracks <- fish_simulation$tracks %>%
  rename(fish_id = path_id, time = datetime, x_true = x, y_true = y) %>%
  select(fish_id, time, x_true, y_true)

# Merge with particle filter estimates
comparison <- pf_results$positions %>%
  left_join(true_tracks, by = c("fish_id", "time")) %>%
  filter(!is.na(x_true)) %>%
  mutate(
    error_m = sqrt((x_mean - x_true)^2 + (y_mean - y_true)^2)
  )

cat("\n=== Position Error Summary ===\n")
cat("Mean error:", round(mean(comparison$error_m, na.rm = TRUE), 1), "m\n")
cat("Median error:", round(median(comparison$error_m, na.rm = TRUE), 1), "m\n")
cat("95th percentile:", round(quantile(comparison$error_m, 0.95, na.rm = TRUE), 1), "m\n")

# Error by fish
error_by_fish <- comparison %>%
  group_by(fish_id) %>%
  summarise(
    mean_error = mean(error_m, na.rm = TRUE),
    median_error = median(error_m, na.rm = TRUE),
    n_steps = n(),
    .groups = "drop"
  )
print(error_by_fish)


# 6. VISUALIZATION =============================================================

raster_df <- as.data.frame(depth_raster, xy = TRUE)

# Plot: estimated vs true positions for each fish
p1 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_path(data = comparison, aes(x = x_true, y = y_true),
            col = "white", alpha = 0.5, linewidth = 0.3) +
  geom_point(data = comparison, aes(x = x_mean, y = y_mean, colour = error_m),
             size = 0.8) +
  scale_colour_viridis_c(option = "magma", name = "Error (m)") +
  facet_wrap(~fish_id) +
  theme_minimal() +
  labs(title = "Particle Filter Position Estimates",
       subtitle = "Points = estimates, white line = true track")

print(p1)

# Plot: position error over time
p2 <- ggplot(comparison, aes(x = time, y = error_m)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(colour = n_detecting), size = 0.8) +
  scale_colour_viridis_c(name = "Detecting\nStations") +
  facet_wrap(~fish_id, scales = "free_x") +
  theme_minimal() +
  labs(title = "Position Error Over Time",
       x = "Time", y = "Error (m)")

print(p2)

# Plot: error vs number of detecting stations
p3 <- ggplot(comparison, aes(x = factor(n_detecting), y = error_m)) +
  geom_boxplot(fill = "steelblue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Position Error by Number of Detecting Stations",
       x = "Detecting Stations", y = "Error (m)")

print(p3)

# Plot: ESS over time
p4 <- ggplot(comparison, aes(x = time, y = ess)) +
  geom_line() +
  geom_hline(yintercept = 1000 * 0.5, linetype = "dashed", colour = "red") +
  geom_point(aes(colour = resampled), size = 0.8) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "red"),
                      name = "Resampled") +
  facet_wrap(~fish_id, scales = "free_x") +
  theme_minimal() +
  labs(title = "Effective Sample Size Over Time",
       subtitle = "Red dashed = resampling threshold",
       x = "Time", y = "ESS")

print(p4)


# 7. BENCHMARK SCALING =========================================================

cat("\n=== Benchmarking particle count scaling ===\n")
for (np in c(100, 500, 1000, 2000)) {
  t <- system.time({
    particle_filter_positioning(
      detection_data = fish_simulation$station_detections,
      station_info = points_regular,
      de_model = logistic_DE$log_model,
      raster = depth_raster,
      n_particles = np,
      step_length_mean = 50, step_length_sd = 30,
      turning_angle_sd = 45, time_step = 180,
      max_distance = 30000,
      fish_id_col = "path_id", time_col = "datetime",
      station_col = "station_id",
      verbose = FALSE
    )
  })
  cat(sprintf("  %4d particles: %.1f s\n", np, t["elapsed"]))
}

cat("\nDone!\n")
