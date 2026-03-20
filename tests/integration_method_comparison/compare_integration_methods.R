# Compare Integration Methods for WADE Positioning
# This script runs the same positioning scenario with all three integration methods
# (subtractive, multiplicative, additive) and produces side-by-side comparison plots.

# 1. SETUP =====================================================================
library(positionR)
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)

# 2. DATA & ARRAY SETUP =======================================================

data("depth_raster")
raster::crs(depth_raster) <- "EPSG:32617"

points_regular <- generate_exact_regular_points(depth_raster, n_points = 40, seed = 123)

station_distances <- calculate_station_distances(
  raster = depth_raster,
  receiver_frame = points_regular,
  max_distance = 30000,
  station_col = "station_id"
)

# 3. DETECTION RANGE MODEL ====================================================

logistic_DE <- create_logistic_curve_depth(
  min_depth = 1, max_depth = 35,
  d50_min_depth = 400, d95_min_depth = 800,
  d50_max_depth = 750, d95_max_depth = 1500,
  plot = FALSE, return_model = TRUE, return_object = TRUE
)

station_distances$DE_pred <- stats::predict(
  logistic_DE$log_model,
  newdata = station_distances %>%
    dplyr::rename(dist_m = cost_distance) %>%
    dplyr::mutate(depth_m = abs(raster_value)),
  type = "response"
)

# 4. FISH SIMULATION ===========================================================

start_time <- as.POSIXct("2025-07-15 12:00:00", tz = "UTC")

fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances,
  n_paths = 5,
  n_steps = 480,
  step_length_mean = 50,
  step_length_sd = 30,
  time_step = 180,
  seed = 123,
  start_time = start_time,
  include_barriers = TRUE
)

# 5. RUN WADE WITH ALL THREE INTEGRATION METHODS ==============================

cat("\n=== Running SUBTRACTIVE integration ===\n")
results_subtractive <- calculate_fish_positions(
  station_detections = fish_simulation$station_detections,
  station_distances_df = station_distances,
  station_info = points_regular,
  de_model = logistic_DE$log_model,
  integration_method = "subtractive",
  non_detection_weight = 0.5,
  max_non_detection_distance = 2000,
  weighting_method = "normalize_stations",
  normalization_method = "min_max",
  fish_id_col = "path_id",
  time_col = "datetime",
  time_aggregation = "day",
  station_col = "station_id",
  include_barriers = TRUE,
  verbose = FALSE
)

cat("\n=== Running MULTIPLICATIVE integration ===\n")
results_multiplicative <- calculate_fish_positions(
  station_detections = fish_simulation$station_detections,
  station_distances_df = station_distances,
  station_info = points_regular,
  de_model = logistic_DE$log_model,
  integration_method = "multiplicative",
  non_detection_weight = 0.5,
  max_non_detection_distance = 2000,
  weighting_method = "normalize_stations",
  normalization_method = "min_max",
  fish_id_col = "path_id",
  time_col = "datetime",
  time_aggregation = "day",
  station_col = "station_id",
  include_barriers = TRUE,
  verbose = FALSE
)

cat("\n=== Running ADDITIVE integration ===\n")
results_additive <- calculate_fish_positions(
  station_detections = fish_simulation$station_detections,
  station_distances_df = station_distances,
  station_info = points_regular,
  de_model = logistic_DE$log_model,
  integration_method = "additive",
  detection_weight = 0.5,
  non_detection_weight = 0.5,
  max_non_detection_distance = 2000,
  weighting_method = "normalize_stations",
  normalization_method = "min_max",
  fish_id_col = "path_id",
  time_col = "datetime",
  time_aggregation = "day",
  station_col = "station_id",
  include_barriers = TRUE,
  verbose = FALSE
)

# 6. COMPARISON PLOTS ==========================================================

fish_select <- 1
time_select <- "2025-07-15"

# --- All 5 panels in a single plot, saved to PNG ---

# Helper to strip redundant per-panel elements
strip_panel <- function(p, title) {
  p + ggtitle(title) +
    theme(
      plot.subtitle = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

p_det <- plot_fish_positions(
  positioning_results = results_subtractive,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "detection",
  prob_threshold = 0.01,
  track_data = fish_simulation$tracks
) %>% strip_panel("Detection Field")

p_nondet <- plot_fish_positions(
  positioning_results = results_subtractive,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "non_detection",
  track_data = fish_simulation$tracks
) %>% strip_panel("Non-Detection Field")

p_sub <- plot_fish_positions(
  positioning_results = results_subtractive,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "integrated",
  prob_threshold = 0.01,
  detection_threshold = 0.01,
  track_data = fish_simulation$tracks
) %>% strip_panel("Integrated: Subtractive")

p_mult <- plot_fish_positions(
  positioning_results = results_multiplicative,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "integrated",
  prob_threshold = 0.01,
  detection_threshold = 0.01,
  track_data = fish_simulation$tracks
) %>% strip_panel("Integrated: Multiplicative")

p_add <- plot_fish_positions(
  positioning_results = results_additive,
  depth_raster_df = depth_raster,
  fish_select = fish_select,
  time_select = time_select,
  plot_type = "integrated",
  prob_threshold = 0.01,
  detection_threshold = 0.01,
  track_data = fish_simulation$tracks
) %>% strip_panel("Integrated: Additive (original)")

# Input fields top row (centred), integrated methods bottom row
combined_plot <- (p_det + p_nondet + plot_spacer()) /
  (p_sub + p_mult + p_add) +
  plot_annotation(
    title = "Integration Method Comparison",
    subtitle = paste("Fish", fish_select, "| Time:", time_select,
                     "| non_detection_weight = 0.5")
  )

combined_plot

# --- Summary statistics ---

cat("\n=== Summary Statistics ===\n\n")

summarise_method <- function(results, method_name) {
  probs <- results$position_probabilities %>%
    filter(fish_id == fish_select)

  # Filter by time_period_label (date string like "2025-07-15")
  if (is.character(time_select) && "time_period_label" %in% names(probs)) {
    probs <- probs %>% filter(grepl(time_select, time_period_label))
  }

  data.frame(
    method = method_name,
    n_cells_above_0 = sum(probs$integrated_prob > 0, na.rm = TRUE),
    n_cells_above_0.1 = sum(probs$integrated_prob > 0.1, na.rm = TRUE),
    n_cells_above_0.5 = sum(probs$integrated_prob > 0.5, na.rm = TRUE),
    mean_prob = round(mean(probs$integrated_prob, na.rm = TRUE), 4),
    max_prob = round(max(probs$integrated_prob, na.rm = TRUE), 4),
    min_prob = round(min(probs$integrated_prob, na.rm = TRUE), 4)
  )
}

summary_df <- bind_rows(
  summarise_method(results_subtractive, "subtractive"),
  summarise_method(results_multiplicative, "multiplicative"),
  summarise_method(results_additive, "additive")
)

print(summary_df)

cat("\nSubtractive and multiplicative should show fewer cells with probability > 0\n")
cat("compared to additive, demonstrating the tighter spatial footprint.\n")


# --- Multi-fish comparison ---

cat("\n=== Multi-fish comparison (all fish, all days) ===\n\n")

multi_summary <- function(results, method_name) {
  probs <- results$position_probabilities
  probs %>%
    group_by(fish_id, time_period) %>%
    summarise(
      n_cells_gt0 = sum(integrated_prob > 0),
      mean_prob = mean(integrated_prob),
      .groups = "drop"
    ) %>%
    summarise(
      method = method_name,
      mean_cells_gt0 = round(mean(n_cells_gt0)),
      sd_cells_gt0 = round(sd(n_cells_gt0)),
      mean_avg_prob = round(mean(mean_prob), 4),
      .groups = "drop"
    )
}

multi_df <- bind_rows(
  multi_summary(results_subtractive, "subtractive"),
  multi_summary(results_multiplicative, "multiplicative"),
  multi_summary(results_additive, "additive")
)

print(multi_df)

# fin.
