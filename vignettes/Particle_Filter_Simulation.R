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
points_regular <- generate_exact_regular_points(depth_raster, n_points = 80, seed = 123)

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
  n_paths = 6,
  n_steps = 480,
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
    return_particles = TRUE,  # Need particles for visualization
    verbose = TRUE
  )
})
cat("Particle filter time:", pf_time["elapsed"], "s\n")
cat("Position estimates:", nrow(pf_results$positions), "\n")
cat("Particle records:", nrow(pf_results$particles), "\n")


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

# Prepare receiver station coordinates for plotting
station_coords <- sf::st_coordinates(points_regular)
stations_df <- data.frame(
  station_id = points_regular$station_id,
  station_x = station_coords[, 1],
  station_y = station_coords[, 2]
)

# Summarize detections per station per fish (for correct faceting)
det_summary_by_fish <- fish_simulation$station_detections %>%
  filter(detected == 1) %>%
  rename(fish_id = path_id) %>%
  group_by(fish_id, station_id) %>%
  summarise(total_dets = n(), .groups = "drop")

# Create per-fish station data: all stations for each fish, with detection counts
fish_ids_plot <- unique(pf_results$positions$fish_id)
stations_plot <- do.call(rbind, lapply(fish_ids_plot, function(fid) {
  s <- stations_df
  s$fish_id <- fid
  s
}))
stations_plot <- stations_plot %>%
  left_join(det_summary_by_fish, by = c("fish_id", "station_id")) %>%
  mutate(total_dets = ifelse(is.na(total_dets), 0, total_dets),
         has_detections = total_dets > 0)


# Compute zoom extent from true tracks + estimates with 15% buffer
all_x <- c(comparison$x_true, comparison$x_mean)
all_y <- c(comparison$y_true, comparison$y_mean)
x_range <- range(all_x, na.rm = TRUE)
y_range <- range(all_y, na.rm = TRUE)
x_buf <- diff(x_range) * 0.15
y_buf <- diff(y_range) * 0.15
zoom_xlim <- c(x_range[1] - x_buf, x_range[2] + x_buf)
zoom_ylim <- c(y_range[1] - y_buf, y_range[2] + y_buf)


# --- Plot 1: Particle cloud + resolved paths + true tracks ---
# Subsample particles to avoid overplotting
# Sample every 10th time step, 50 particles each
pf_times <- unique(pf_results$particles$time)
sampled_times <- pf_times[seq(1, length(pf_times), by = 10)]
particles_sub <- pf_results$particles %>%
  filter(time %in% sampled_times) %>%
  group_by(fish_id, time) %>%
  slice_sample(n = 500) %>%
  ungroup()

p1 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  # Particle cloud (light, transparent)
  geom_point(data = particles_sub, aes(x = x, y = y),
             colour = "yellow", alpha = 0.1, size = 0.3) +
  # Estimated track (white)
  geom_path(data = comparison, aes(x = x_mean, y = y_mean),
            colour = "white", linewidth = 0.5, alpha = 0.8) +
  geom_point(data = comparison, aes(x = x_mean, y = y_mean),
             colour = "white", size = 0.5, alpha = 0.6) +
  # True track (green)
  geom_path(data = comparison, aes(x = x_true, y = y_true),
            colour = "green3", linewidth = 0.4, alpha = 0.7) +
  # Receiver stations
  geom_point(data = stations_plot %>% filter(!has_detections),
             aes(x = station_x, y = station_y),
             colour = "red", size = 1, shape = 4) +
  geom_point(data = stations_plot %>% filter(has_detections),
             aes(x = station_x, y = station_y, size = total_dets),
             colour = "orange", fill = NA, shape = 21, stroke = 0.8) +
  scale_size_continuous(range = c(1, 5), name = "Detections") +
  facet_wrap(~fish_id, ncol = 2) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = "Particle Filter: Particle Cloud + Resolved Paths",
       subtitle = "Yellow = particles | White = estimated path | Green = true path")

print(p1)


# --- Plot 2: Estimated vs true with error colouring ---
p2 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  # Estimated positions coloured by error
  geom_point(data = comparison, aes(x = x_mean, y = y_mean, colour = error_m),
             size = 1) +
  scale_colour_viridis_c(option = "magma", name = "Error (m)") +
  # Uncertainty bars (95% CI ~ 2*SD)
  geom_segment(data = comparison,
               aes(x = x_mean - 2*x_sd, xend = x_mean + 2*x_sd,
                   y = y_mean, yend = y_mean),
               colour = "red", alpha = 0.3, linewidth = 0.2) +
  geom_segment(data = comparison,
               aes(x = x_mean, xend = x_mean,
                   y = y_mean - 2*y_sd, yend = y_mean + 2*y_sd),
               colour = "red", alpha = 0.3, linewidth = 0.2) +
  # True track
  geom_path(data = comparison, aes(x = x_true, y = y_true),
            colour = "green3", linewidth = 0.3, alpha = 0.6) +
  # Stations
  geom_point(data = stations_plot, aes(x = station_x, y = station_y),
             colour = "grey80", size = 0.5, shape = 3) +
  facet_wrap(~fish_id, ncol = 2) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = "Position Estimates with Error and Uncertainty",
       subtitle = "Green = true track | Red bars = 95% CI | Colour = position error")

print(p2)


# --- Plot 3: Particle cloud at selected time steps ---
# Show particle spread at a few key time steps for one fish
fish1_times <- comparison %>%
  filter(fish_id == comparison$fish_id[1]) %>%
  arrange(time)
# Select ~6 evenly spaced time steps
time_idx <- round(seq(1, nrow(fish1_times), length.out = 6))
selected_times <- fish1_times$time[time_idx]

particles_snapshot <- pf_results$particles %>%
  filter(fish_id == comparison$fish_id[1], time %in% selected_times)
positions_snapshot <- comparison %>%
  filter(fish_id == comparison$fish_id[1], time %in% selected_times)

p3 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  # Particles coloured by weight
  geom_point(data = particles_snapshot,
             aes(x = x, y = y, colour = weight),
             size = 0.5, alpha = 0.4) +
  scale_colour_viridis_c(option = "magma", name = "Weight") +
  # Weighted mean position (white)
  geom_point(data = positions_snapshot,
             aes(x = x_mean, y = y_mean),
             colour = "white", size = 2, shape = 4, stroke = 1.5) +
  # True position (green)
  geom_point(data = positions_snapshot,
             aes(x = x_true, y = y_true),
             colour = "green3", size = 2, shape = 3, stroke = 1.5) +
  # Stations
  geom_point(data = stations_plot, aes(x = station_x, y = station_y),
             colour = "grey80", size = 0.5, shape = 3) +
  facet_wrap(~time, ncol = 2) +
  coord_sf(xlim = range(c(particles_snapshot$x, positions_snapshot$x_mean), na.rm = TRUE) + c(-500, 500),
           ylim = range(c(particles_snapshot$y, positions_snapshot$y_mean), na.rm = TRUE) + c(-500, 500)) +
  theme_minimal() +
  labs(title = paste("Particle Distributions at Selected Time Steps - Fish",
                     comparison$fish_id[1]),
       subtitle = "White cross = estimate | Green cross = true position | Colour = particle weight")

print(p3)


# --- Plot 4: Position error over time ---
p4 <- ggplot(comparison, aes(x = time, y = error_m)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(colour = n_detecting), size = 0.8) +
  scale_colour_viridis_c(name = "Detecting\nStations") +
  facet_wrap(~fish_id, scales = "free_x") +
  theme_minimal() +
  labs(title = "Position Error Over Time",
       x = "Time", y = "Error (m)")

print(p4)


# --- Plot 5: ESS over time ---
p5 <- ggplot(comparison, aes(x = time, y = ess)) +
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

print(p5)


# --- Plot 6: Error by detecting stations ---
p6 <- ggplot(comparison, aes(x = factor(n_detecting), y = error_m)) +
  geom_boxplot(fill = "steelblue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Position Error by Number of Detecting Stations",
       x = "Detecting Stations", y = "Error (m)")

print(p6)


# 7. TWO-FILTER SMOOTHER ======================================================

# The smoother runs a backward filter and combines with the forward pass,
# incorporating both past and future observations at each time step.

cat("\n=== Running Two-Filter Smoother ===\n")
smooth_time <- system.time({
  smooth_results <- pf_smooth(
    pf_results = pf_results,
    detection_data = fish_simulation$station_detections,
    station_info = points_regular,
    de_model = logistic_DE$log_model,
    raster = depth_raster,
    step_length_mean = 50,
    step_length_sd = 30,
    turning_angle_sd = 45,
    time_step = 180,
    max_distance = 30000,
    fish_id_col = "path_id",
    time_col = "datetime",
    station_col = "station_id",
    return_particles = TRUE,
    verbose = TRUE
  )
})
cat("Smoother time:", smooth_time["elapsed"], "s\n")

# Compare forward vs smoothed error
smooth_comparison <- smooth_results$positions %>%
  left_join(true_tracks, by = c("fish_id", "time")) %>%
  filter(!is.na(x_true)) %>%
  mutate(error_m = sqrt((x_mean - x_true)^2 + (y_mean - y_true)^2))

cat("\n=== Forward vs Smoothed Error ===\n")
cat("Forward  - mean:", round(mean(comparison$error_m, na.rm = TRUE), 1),
    "m, median:", round(median(comparison$error_m, na.rm = TRUE), 1), "m\n")
cat("Smoothed - mean:", round(mean(smooth_comparison$error_m, na.rm = TRUE), 1),
    "m, median:", round(median(smooth_comparison$error_m, na.rm = TRUE), 1), "m\n")
cat("Improvement:", round((1 - mean(smooth_comparison$error_m, na.rm = TRUE) /
    mean(comparison$error_m, na.rm = TRUE)) * 100, 1), "%\n")

# Plot: forward vs smoothed positions
p7 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  # True track
  geom_path(data = comparison, aes(x = x_true, y = y_true),
            colour = "green3", linewidth = 0.4, alpha = 0.7) +
  # Smoothed estimate (white)
  geom_path(data = smooth_comparison, aes(x = x_mean, y = y_mean),
            colour = "white", linewidth = 0.5, alpha = 0.8) +
  # Forward estimate (red)
  geom_path(data = comparison, aes(x = x_mean, y = y_mean),
            colour = "red", linewidth = 0.3, alpha = 0.5) +
  # Stations
  geom_point(data = stations_plot %>% filter(!has_detections),
             aes(x = station_x, y = station_y),
             colour = "red", size = 1, shape = 4) +
  geom_point(data = stations_plot %>% filter(has_detections),
             aes(x = station_x, y = station_y, size = total_dets),
             colour = "orange", fill = NA, shape = 21, stroke = 0.8) +
  scale_size_continuous(range = c(1, 5), name = "Detections") +
  facet_wrap(~fish_id, ncol = 2) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = "Forward vs Smoothed Position Estimates",
       subtitle = "Green = true | Red = forward | White = smoothed")

print(p7)


# 8. UTILIZATION DISTRIBUTION =================================================

cat("\n=== Computing Utilization Distributions ===\n")
ud_results <- pf_utilization_distribution(
  pf_results = smooth_results,
  raster = depth_raster,
  contour_levels = c(0.5, 0.95),
  by_time = FALSE,
  verbose = TRUE
)

cat("\nHome Range Estimates:\n")
print(ud_results$home_ranges)

# Plot UD for first fish
first_fish_id <- as.character(unique(smooth_results$positions$fish_id)[1])
ud_raster_1 <- ud_results$ud_rasters[[first_fish_id]]
ud_df <- as.data.frame(ud_raster_1, xy = TRUE)
names(ud_df)[3] <- "prob"

p8 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = ud_df %>% filter(!is.na(prob) & prob > 0),
              aes(x = x, y = y, fill = prob), alpha = 0.8) +
  scale_fill_viridis_c(option = "magma", name = "P(use)") +
  geom_path(data = comparison %>% filter(fish_id == as.integer(first_fish_id)),
            aes(x = x_true, y = y_true),
            colour = "green3", linewidth = 0.4, alpha = 0.7) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = paste("Utilization Distribution - Fish", first_fish_id),
       subtitle = "Green = true track | Colour = probability of use")

print(p8)

# Plot 50% and 95% UD contours for all fish
# Classify UD cells into contour bands
ud_contour_list <- list()
for (fid in names(ud_results$ud_rasters)) {
  ud_r <- ud_results$ud_rasters[[fid]]
  ud_vals <- as.data.frame(ud_r, xy = TRUE)
  names(ud_vals)[3] <- "prob"
  ud_vals <- ud_vals %>% filter(!is.na(prob) & prob > 0)

  # Sort by probability descending, compute cumulative sum
  ud_vals <- ud_vals %>% arrange(desc(prob)) %>%
    mutate(cum_prob = cumsum(prob),
           contour = case_when(
             cum_prob <= 0.50 ~ "Core (50%)",
             cum_prob <= 0.95 ~ "Home range (95%)",
             TRUE ~ NA_character_
           ),
           fish_id = as.integer(fid))
  ud_contour_list[[fid]] <- ud_vals %>% filter(!is.na(contour))
}
ud_contours_df <- do.call(rbind, ud_contour_list)

p9 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = ud_contours_df,
              aes(x = x, y = y, fill = contour), alpha = 0.7) +
  scale_fill_manual(values = c("Core (50%)" = "red", "Home range (95%)" = "orange"),
                    name = "UD Contour") +
  # True track
  geom_path(data = comparison, aes(x = x_true, y = y_true),
            colour = "green3", linewidth = 0.3, alpha = 0.7) +
  # Smoothed track
  geom_path(data = smooth_comparison, aes(x = x_mean, y = y_mean),
            colour = "white", linewidth = 0.3, alpha = 0.5) +
  facet_wrap(~fish_id, ncol = 2) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = "Utilization Distribution Contours (50% and 95%)",
       subtitle = "Green = true track | White = smoothed estimate")

print(p9)


# 9. DETECTION GAP ANALYSIS ====================================================

cat("\n=== Detection Gap Analysis ===\n")
gap_results <- pf_analyze_gaps(
  pf_results = smooth_results,
  true_tracks = fish_simulation$tracks %>%
    rename(fish_id = path_id, time = datetime, x_true = x, y_true = y) %>%
    select(fish_id, time, x_true, y_true)
)

cat("Total gaps identified:", nrow(gap_results$gaps), "\n")
if (nrow(gap_results$gaps) > 0) {
  cat("Gap duration range:",
      round(min(gap_results$gaps$gap_duration_min), 1), "-",
      round(max(gap_results$gaps$gap_duration_min), 1), "min\n")
  cat("Mean max error during gaps:",
      round(mean(gap_results$gaps$max_error_in_gap, na.rm = TRUE), 1), "m\n")
}

# Plot: error vs time since last detection
ts <- gap_results$time_series
p10 <- ggplot(ts, aes(x = time_since_last_detection_sec / 60, y = error_m)) +
  geom_point(alpha = 0.3, size = 0.5, colour = "steelblue") +
  geom_smooth(method = "loess", colour = "red", se = TRUE, linewidth = 0.8) +
  theme_minimal() +
  labs(title = "Position Error vs Time Since Last Detection",
       x = "Time since last detection (min)", y = "Error (m)")

print(p10)

# Plot: error trajectories per gap, coloured by gap duration
if (nrow(gap_results$gaps) > 0) {
  gap_ts <- ts %>% filter(!is.na(gap_id))
  gap_ts <- gap_ts %>%
    left_join(gap_results$gaps %>% select(gap_id, gap_duration_min),
              by = "gap_id")

  p11 <- ggplot(gap_ts, aes(x = time_since_last_detection_sec / 60,
                              y = error_m,
                              group = gap_id,
                              colour = gap_duration_min)) +
    geom_line(alpha = 0.5, linewidth = 0.4) +
    scale_colour_viridis_c(name = "Gap\nduration\n(min)") +
    theme_minimal() +
    labs(title = "Error Growth During Detection Gaps",
         subtitle = "Each line = one gap event",
         x = "Time since last detection (min)", y = "Error (m)")

  print(p11)
}


# 10. BENCHMARK SCALING =========================================================

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
