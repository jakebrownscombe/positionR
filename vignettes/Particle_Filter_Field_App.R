# Particle Filter Positioning - Field Data Application
# Jake Brownscombe
# 2025

# This script applies particle filter positioning to real acoustic telemetry data
# (walleye in Stoney Lake). It uses the same data as the WADE field application
# but estimates continuous movement paths via Sequential Monte Carlo.

# 1. SETUP AND LIBRARIES =====================================================
library(positionR)
library(raster)
library(sf)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(patchwork)

set.seed(123)


# 2. DATA LOADING ============================================================

data("stoney_fish_detections")
data("stoney_rx_deploy")
data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617"


# 3. SPATIAL SETUP - STATION DISTANCES =======================================

station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = stoney_rx_deploy,
  max_distance = 30000,
  station_col = "station_id"
)


# 4. DETECTION RANGE MODEL ===================================================

logistic_DE <- create_logistic_curve_depth(
  min_depth = 1,
  max_depth = 35,
  d50_min_depth = 400,
  d95_min_depth = 800,
  d50_max_depth = 750,
  d95_max_depth = 1500,
  plot = TRUE,
  return_model = TRUE,
  return_object = TRUE
)

station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = station_distances %>%
    dplyr::rename(dist_m = cost_distance) %>%
    dplyr::mutate(depth_m = abs(raster_value)),
  type = "response"
)


# 5. PREPARE DETECTION DATA FOR PARTICLE FILTER ==============================

# Select fish and time period (same as WADE field app)
selected_fish <- "Walleye-1504321"
start_time <- as.POSIXct("2023-10-15 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2023-10-29 23:00:00", tz = "UTC")

# Prepare data at exact detection timestamps (no temporal aggregation)
pf_data <- prepare_detection_data_for_pf(
  fish_detections = stoney_fish_detections,
  station_deployments = stoney_rx_deploy,
  selected_fish_id = selected_fish,
  start_time = start_time,
  end_time = end_time,
  min_gap_sec = 60  # Merge detections within 60s
)

head(pf_data)
cat("Unique time steps:", length(unique(pf_data$time)), "\n")


# 6. PARTICLE FILTER POSITIONING =============================================

cat("\n=== Running Particle Filter ===\n")
pf_time <- system.time({
  pf_results <- particle_filter_positioning(
    detection_data = pf_data,
    station_info = stoney_rx_deploy,
    de_model = logistic_DE$log_model,
    raster = depth_raster,
    n_particles = 1000,
    step_length_mean = 50,
    step_length_sd = 30,
    turning_angle_sd = 45,
    time_step = 180,
    max_distance = 30000,
    fish_id_col = "fish_id",
    time_col = "time",
    station_col = "station_id",
    return_particles = TRUE,
    verbose = TRUE
  )
})
cat("Particle filter time:", pf_time["elapsed"], "s\n")
cat("Position estimates:", nrow(pf_results$positions), "\n")


# 7. TWO-FILTER SMOOTHER ====================================================

cat("\n=== Running Smoother ===\n")
smooth_time <- system.time({
  smooth_results <- pf_smooth(
    pf_results = pf_results,
    detection_data = pf_data,
    station_info = stoney_rx_deploy,
    de_model = logistic_DE$log_model,
    raster = depth_raster,
    step_length_mean = 50,
    step_length_sd = 30,
    turning_angle_sd = 45,
    time_step = 180,
    max_distance = 30000,
    fish_id_col = "fish_id",
    time_col = "time",
    station_col = "station_id",
    return_particles = TRUE,
    verbose = TRUE
  )
})
cat("Smoother time:", smooth_time["elapsed"], "s\n")


# 8. VISUALIZATION ==========================================================

raster_df <- as.data.frame(depth_raster, xy = TRUE)

# Station coordinates for plotting
station_coords <- sf::st_coordinates(stoney_rx_deploy)
stations_df <- data.frame(
  station_id = stoney_rx_deploy$station_id,
  station_x = station_coords[, 1],
  station_y = station_coords[, 2]
)

# Detection summary for station sizing
det_summary <- pf_data %>%
  filter(detected == 1) %>%
  group_by(station_id) %>%
  summarise(total_dets = n(), .groups = "drop")
stations_plot <- stations_df %>%
  left_join(det_summary, by = "station_id") %>%
  mutate(total_dets = ifelse(is.na(total_dets), 0, total_dets),
         has_detections = total_dets > 0)

# Zoom extent from estimated positions
pos <- smooth_results$positions
x_range <- range(pos$x_mean, na.rm = TRUE)
y_range <- range(pos$y_mean, na.rm = TRUE)
x_buf <- diff(x_range) * 0.15; y_buf <- diff(y_range) * 0.15
zoom_xlim <- c(x_range[1] - x_buf, x_range[2] + x_buf)
zoom_ylim <- c(y_range[1] - y_buf, y_range[2] + y_buf)

# Plot: smoothed path with receiver array
p1 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  # Smoothed path
  geom_path(data = pos, aes(x = x_mean, y = y_mean),
            colour = "white", linewidth = 0.5, alpha = 0.8) +
  geom_point(data = pos, aes(x = x_mean, y = y_mean),
             colour = "white", size = 0.3, alpha = 0.4) +
  # Stations sized by total detections
  geom_point(data = stations_plot,
             aes(x = station_x, y = station_y, size = total_dets + 1),
             colour = "orange", fill = NA, shape = 21, stroke = 0.8) +
  scale_size_continuous(range = c(0.5, 5), name = "Detections",
                        breaks = c(1, 50, 200, 500)) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = paste("Particle Filter Path Estimate -", selected_fish),
       subtitle = "White = smoothed path | Orange circles = receivers (sized by detections)")

print(p1)


# Plot: path coloured by uncertainty
p2 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  geom_point(data = pos, aes(x = x_mean, y = y_mean, colour = x_sd + y_sd),
             size = 0.8) +
  scale_colour_viridis_c(option = "magma", name = "Uncertainty\n(SD, m)") +
  geom_point(data = stations_plot, aes(x = station_x, y = station_y),
             colour = "grey80", size = 0.5, shape = 3) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = "Position Uncertainty",
       subtitle = "Colour = combined x + y standard deviation")

print(p2)


# 9. UTILIZATION DISTRIBUTIONS ==============================================

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

# Plot UD contours
fish_id_str <- as.character(selected_fish)
ud_r <- ud_results$ud_rasters[[fish_id_str]]
ud_vals <- as.data.frame(ud_r, xy = TRUE)
names(ud_vals)[3] <- "prob"
ud_vals <- ud_vals %>% filter(!is.na(prob) & prob > 0) %>%
  arrange(desc(prob)) %>%
  mutate(cum_prob = cumsum(prob),
         contour = case_when(
           cum_prob <= 0.50 ~ "Core (50%)",
           cum_prob <= 0.95 ~ "Home range (95%)",
           TRUE ~ NA_character_))
ud_contours <- ud_vals %>% filter(!is.na(contour))

p3 <- ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradient(low = "blue4", high = "cornflowerblue",
                      na.value = "transparent", name = "Depth (m)") +
  new_scale_fill() +
  geom_raster(data = ud_contours,
              aes(x = x, y = y, fill = contour), alpha = 0.7) +
  scale_fill_manual(values = c("Core (50%)" = "red", "Home range (95%)" = "orange"),
                    name = "UD Contour") +
  geom_path(data = pos, aes(x = x_mean, y = y_mean),
            colour = "white", linewidth = 0.3, alpha = 0.5) +
  geom_point(data = stations_plot %>% filter(has_detections),
             aes(x = station_x, y = station_y),
             colour = "yellow", size = 1) +
  coord_sf(xlim = zoom_xlim, ylim = zoom_ylim) +
  theme_minimal() +
  labs(title = paste("Utilization Distribution -", selected_fish),
       subtitle = "White = smoothed path | Yellow = detecting stations")

print(p3)


# 10. DETECTION GAP ANALYSIS ================================================

cat("\n=== Detection Gap Analysis ===\n")
gap_results <- pf_analyze_gaps(
  pf_results = smooth_results,
  true_tracks = NULL  # No ground truth for field data
)

cat("Total gaps:", nrow(gap_results$gaps), "\n")
if (nrow(gap_results$gaps) > 0) {
  cat("Gap duration range:",
      round(min(gap_results$gaps$gap_duration_min), 1), "-",
      round(max(gap_results$gaps$gap_duration_min), 1), "min\n")
  cat("Mean uncertainty during gaps:",
      round(mean(gap_results$gaps$mean_uncertainty, na.rm = TRUE), 1), "m\n")
}

# Plot: uncertainty vs time since last detection
ts <- gap_results$time_series
p4 <- ggplot(ts, aes(x = time_since_last_detection_sec / 60,
                       y = x_sd + y_sd)) +
  geom_point(alpha = 0.3, size = 0.5, colour = "steelblue") +
  geom_smooth(method = "loess", colour = "red", se = TRUE, linewidth = 0.8) +
  theme_minimal() +
  labs(title = "Position Uncertainty vs Time Since Last Detection",
       x = "Time since last detection (min)",
       y = "Uncertainty (x_sd + y_sd, m)")

print(p4)

# Plot: ESS over time
p5 <- ggplot(pos, aes(x = time, y = ess)) +
  geom_line(alpha = 0.5) +
  geom_point(aes(colour = n_detecting > 0), size = 0.5) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "green3"),
                      name = "Detected", labels = c("No", "Yes")) +
  geom_hline(yintercept = 1000 * 0.5, linetype = "dashed", colour = "red") +
  theme_minimal() +
  labs(title = "Effective Sample Size Over Time",
       subtitle = "Red dashed = resampling threshold",
       x = "Time", y = "ESS")

print(p5)

cat("\nDone!\n")
