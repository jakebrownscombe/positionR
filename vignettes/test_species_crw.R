# ============================================================================
# Test: Species-Specific CRW and Particle Filter
# Author: Jacob Brownscombe
# Description: Tests mode = "species_empirical" for empirical resampling of
#   step lengths and turn angles in simulate_fish_tracks() and
#   particle_filter_positioning().
# ============================================================================

# ---- Libraries ----

library(positionR)
library(raster)
library(sf)
library(ggplot2)
library(patchwork)
library(data.table)

set.seed(123)

# ---- Setup (same as existing vignettes) ----

cat("=== Setup ===\n")

data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617"

# Receiver array
points_regular <- generate_exact_regular_points(depth_raster, n_points = 80, seed = 123)

# Station distances
station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = points_regular,
  max_distance = 30000,
  station_col = "station_id"
)

# DE model
logistic_DE <- create_logistic_curve_depth(
  min_depth = 1, max_depth = 35,
  d50_min_depth = 400, d95_min_depth = 800,
  d50_max_depth = 750, d95_max_depth = 1500,
  plot = FALSE, return_model = TRUE, return_object = TRUE
)

station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = data.frame(
    dist_m = station_distances$cost_distance,
    depth_m = abs(station_distances$raster_value)
  ),
  type = "response"
)

# ---- 1. Species-Specific CRW Simulations ----

cat("\n=== Simulating species-specific tracks ===\n")

species_list <- c("Walleye", "Smallmouth Bass", "Muskellunge", "Generic")
start_time <- as.POSIXct("2025-07-15 12:00:00", tz = "UTC")

sim_results <- list()
for (sp in species_list) {
  cat("\n--- ", sp, " ---\n")
  sim_results[[sp]] <- simulate_fish_tracks(
    raster = depth_raster,
    station_distances = station_distances,
    n_paths = 3,
    n_steps = 480,
    mode = "species_empirical", species = sp,
    seed = 42,
    start_time = start_time,
    include_barriers = TRUE
  )
}

# ---- 2. Validate Step Length and Turn Angle Distributions ----

cat("\n=== Validating distributions ===\n")

data("crw_empirical_distributions")

all_tracks <- rbindlist(lapply(names(sim_results), function(sp) {
  dt <- as.data.table(sim_results[[sp]]$tracks)
  dt[, species := sp]
  dt
}))

# Compute step lengths and turn angles from simulated tracks
setorder(all_tracks, species, path_id, step)
all_tracks[, `:=`(
  dx = x - data.table::shift(x, 1),
  dy = y - data.table::shift(y, 1)
), by = .(species, path_id)]
all_tracks[, dist_m := sqrt(dx^2 + dy^2)]
all_tracks[, bearing := atan2(dx, dy)]
all_tracks[, turn_angle := bearing - data.table::shift(bearing, 1), by = .(species, path_id)]
all_tracks[, turn_angle := atan2(sin(turn_angle), cos(turn_angle))]

sim_moves <- all_tracks[!is.na(dist_m) & dist_m > 0]

# Build comparison data
plot_list_step <- list()
plot_list_ta <- list()

for (sp in species_list) {
  emp <- crw_empirical_distributions[[sp]]

  emp_dt <- data.table(value = emp$step_lengths, source = "Empirical (stored)")
  sim_dt <- data.table(value = sim_moves[species == sp]$dist_m, source = "Simulated")
  step_dt <- rbind(emp_dt, sim_dt)
  step_cap <- quantile(emp$step_lengths, 0.99)

  plot_list_step[[sp]] <- ggplot(step_dt[value < step_cap], aes(x = value, fill = source)) +
    geom_density(alpha = 0.5, colour = NA) +
    scale_fill_manual(values = c("Empirical (stored)" = "steelblue", "Simulated" = "coral")) +
    labs(title = sp, x = "Step length (m)", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")

  emp_ta <- data.table(value = emp$turn_angles * 180 / pi, source = "Empirical (stored)")
  sim_ta <- data.table(
    value = sim_moves[species == sp & !is.na(turn_angle)]$turn_angle * 180 / pi,
    source = "Simulated"
  )
  ta_dt <- rbind(emp_ta, sim_ta)

  plot_list_ta[[sp]] <- ggplot(ta_dt, aes(x = value, fill = source)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 10,
                    alpha = 0.5, position = "identity") +
    scale_fill_manual(values = c("Empirical (stored)" = "steelblue", "Simulated" = "coral")) +
    labs(title = sp, x = "Turn angle (degrees)", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")
}

p_steps <- wrap_plots(plot_list_step, ncol = 2) +
  plot_annotation(title = "Step Length: Empirical vs Simulated")
p_turns <- wrap_plots(plot_list_ta, ncol = 2) +
  plot_annotation(title = "Turn Angle: Empirical vs Simulated")
p_validation <- p_steps / p_turns

p_validation

# ---- 3. Simulated Track Maps ----

raster_df <- as.data.frame(depth_raster, xy = TRUE)
station_coords <- sf::st_coordinates(points_regular)
stations_df <- data.frame(
  station_x = station_coords[, 1],
  station_y = station_coords[, 2]
)

track_map_list <- list()
for (sp in species_list) {
  sp_tracks <- as.data.table(sim_results[[sp]]$tracks)

  track_map_list[[sp]] <- ggplot() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                        na.value = "transparent", guide = "none") +
    geom_point(data = stations_df, aes(x = station_x, y = station_y),
               colour = "grey80", size = 0.5, shape = 3) +
    geom_path(data = sp_tracks, aes(x = x, y = y, colour = factor(path_id)),
              linewidth = 0.4, alpha = 0.7) +
    scale_colour_brewer(palette = "Set1", guide = "none") +
    coord_equal() +
    labs(title = sp) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text = element_blank())
}

p_track_maps <- wrap_plots(track_map_list, ncol = 2) +
  plot_annotation(title = "Simulated Tracks by Species (Empirical CRW)")

p_track_maps

# ---- 4. Step Length and Turn Angle Distribution Plots ----

dist_step_list <- list()
dist_ta_list <- list()

for (sp in species_list) {
  sp_moves <- sim_moves[species == sp]

  dist_step_list[[sp]] <- ggplot(sp_moves[dist_m < quantile(dist_m, 0.99)],
                                  aes(x = dist_m)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                    fill = "steelblue", colour = "white", alpha = 0.7) +
    geom_density(colour = "darkblue", linewidth = 0.8) +
    labs(title = sp, x = "Step length (m)", y = "Density") +
    theme_minimal()

  ta_data <- sp_moves[!is.na(turn_angle)]
  dist_ta_list[[sp]] <- ggplot(ta_data, aes(x = turn_angle * 180 / pi)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 10,
                    fill = "coral", colour = "white", alpha = 0.7) +
    labs(title = sp, x = "Turn angle (degrees)", y = "Density") +
    theme_minimal()
}

p_dist <- (wrap_plots(dist_step_list, ncol = 2) +
             plot_annotation(title = "Simulated Step Length Distributions")) /
  (wrap_plots(dist_ta_list, ncol = 2) +
     plot_annotation(title = "Simulated Turn Angle Distributions"))

p_dist

# ---- 5. Species-Specific Particle Filter ----

cat("\n=== Running particle filter by species ===\n")

pf_results_list <- list()
comparison_list <- list()

for (sp in c("Walleye", "Smallmouth Bass", "Muskellunge")) {
  cat("\n--- Particle filter:", sp, "---\n")

  sim <- sim_results[[sp]]

  pf_result <- particle_filter_positioning(
    detection_data = sim$station_detections,
    station_info = points_regular,
    de_model = logistic_DE$log_model,
    raster = depth_raster,
    n_particles = 500,
    mode = "species_empirical", species = sp,
    max_distance = 30000,
    fish_id_col = "path_id",
    time_col = "datetime",
    station_col = "station_id",
    position_method = "mean",
    return_particles = TRUE,
    verbose = TRUE
  )

  pf_results_list[[sp]] <- pf_result

  # Compare with ground truth
  true_tracks <- as.data.table(sim$tracks)
  setnames(true_tracks, c("x", "y"), c("x_true", "y_true"))

  comp <- merge(
    as.data.table(pf_result$positions),
    true_tracks[, .(fish_id = path_id, time = datetime, x_true, y_true)],
    by = c("fish_id", "time"), all.x = TRUE
  )
  comp <- comp[!is.na(x_true)]
  comp[, error_m := sqrt((x_mean - x_true)^2 + (y_mean - y_true)^2)]
  comp[, species := sp]

  comparison_list[[sp]] <- comp

  cat("  Mean error:", round(mean(comp$error_m, na.rm = TRUE), 1), "m\n")
  cat("  Median error:", round(median(comp$error_m, na.rm = TRUE), 1), "m\n")
  cat("  95th percentile:", round(quantile(comp$error_m, 0.95, na.rm = TRUE), 1), "m\n")
}

# ---- 6. Particle Filter Track Maps (True vs Fitted) ----

cat("\n=== Creating particle filter maps ===\n")

pf_map_list <- list()
for (sp in names(comparison_list)) {
  comp <- comparison_list[[sp]]

  # Subsample particles for plotting (first fish only)
  first_fish <- comp$fish_id[1]
  particles_dt <- as.data.table(pf_results_list[[sp]]$particles)
  particles_sub <- particles_dt[fish_id == first_fish]
  particles_sub <- particles_sub[sample(.N, min(.N, 20000))]

  comp_fish <- comp[fish_id == first_fish]

  pf_map_list[[sp]] <- ggplot() +
    geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
    scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                        na.value = "transparent", guide = "none") +
    geom_point(data = particles_sub, aes(x = x, y = y),
               colour = "yellow", alpha = 0.05, size = 0.2) +
    geom_path(data = comp_fish, aes(x = x_true, y = y_true),
              colour = "green3", linewidth = 0.5, alpha = 0.8) +
    geom_path(data = comp_fish, aes(x = x_mean, y = y_mean),
              colour = "white", linewidth = 0.5, alpha = 0.8) +
    geom_point(data = stations_df, aes(x = station_x, y = station_y),
               colour = "grey80", size = 0.5, shape = 3) +
    coord_equal() +
    labs(title = paste0(sp, " (green = true, white = fitted)")) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.text = element_blank())
}

p_pf_maps <- wrap_plots(pf_map_list, ncol = 1) +
  plot_annotation(title = "Particle Filter Positioning: True vs Fitted Tracks")

p_pf_maps

# ---- 7. Species Comparison ----

cat("\n=== Species Comparison ===\n")

all_comp <- rbindlist(comparison_list)

summary_dt <- all_comp[, .(
  mean_error = round(mean(error_m, na.rm = TRUE), 1),
  median_error = round(median(error_m, na.rm = TRUE), 1),
  p95_error = round(quantile(error_m, 0.95, na.rm = TRUE), 1),
  n_positions = .N
), by = species]

print(summary_dt)

p_error <- ggplot(all_comp, aes(x = species, y = error_m, fill = species)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  coord_cartesian(ylim = c(0, quantile(all_comp$error_m, 0.95))) +
  labs(title = "Particle Filter Positioning Error by Species (Empirical CRW)",
       x = "Species", y = "Error (m)") +
  theme_minimal() +
  theme(legend.position = "none")

p_error

cat("\n=== Done ===\n")
