# positionR 1.6.0

## WADE default changes

* **`integration_method` default changed from `"subtractive"` to `"additive"`**
  in `calculate_fish_positions()`. Empirical per-track comparisons showed that
  subtractive integration combined with per-group [0, 1] rescaling (introduced
  in v1.4.0) over-penalises legitimate detection zones and produces visually
  and metrically incorrect WADE fields. The additive formulation
  `(det * detection_weight) + ((1 - nondet) * non_detection_weight)` restores
  graded, detection-centred probability. The subtractive method remains
  available by setting `integration_method = "subtractive"` explicitly.

* **New parameter `min_detection_prob` (default 0.05)** in
  `calculate_fish_positions()` drops cells with
  `weighted_mean_DE_normalized <= min_detection_prob` from
  `position_probabilities` before returning. Under the additive formula,
  cells far from any receiver would otherwise receive positive
  `integrated_prob` via the `(1 - nondet) * non_detection_weight` term alone,
  inflating downstream point clouds, UDs, and habitat-selection estimates.
  Set `min_detection_prob = 0` to disable.

* **`non_detection_weight` default changed from `1` to `0.5`** so the default
  additive integration satisfies the required sum-to-1 constraint
  (`detection_weight + non_detection_weight == 1`). Workflows that explicitly
  set `integration_method = "subtractive"` and pass `non_detection_weight = 1`
  are unaffected. Under subtractive/multiplicative methods `detection_weight`
  is ignored, so this default change only affects additive integration.

## Migration

To reproduce pre-1.6.0 behaviour exactly, pass
`integration_method = "subtractive", non_detection_weight = 1,
min_detection_prob = 0`. Workflows that already set these parameters
explicitly are unaffected by the default changes.


# positionR 1.5.0

## New Features

* **Effort-balanced non-detection weighting**: New `scale_non_detections`
  parameter (default TRUE) in `calculate_fish_positions()` scales non-detection
  evidence by total detections per fish/time period, balancing the
  detection-weighted mean that favours high-count stations.

## Performance

* **data.table conversion**: All internal WADE helper functions converted from
  dplyr to data.table for ~10x speedup on large datasets. Affects
  `aggregate_detections_for_prediction()`, `create_non_detections()`,
  `aggregate_non_detections()`, `aggregate_probability()`,
  `normalize_DE_by_station()`, and information weighting helpers.

* **Vectorized barrier detection**: New `check_line_crosses_barrier_vectorized()`
  calls `raster::extract()` once per station instead of once per cell (~20-100x
  faster for barrier checking in `calculate_station_distances()`).

* **Optimized reshaping**: `tidyr::pivot_longer()` replaced with
  `data.table::melt()` in `calculate_station_distances()`.

## Improvements

* Reduced console output verbosity: `calculate_station_distances()` progress
  line only in interactive sessions; removed per-station detection rate printout
  from `analyze_detection_performance()`.

* New dependency: `data.table` added to Imports.

---

# positionR 1.4.0

## Breaking Changes

* **New default integration method**: `calculate_fish_positions()` now defaults to
  `integration_method = "subtractive"`, which produces tighter position estimates
  anchored to actual detection zones. The previous additive method inflated spatial
  footprint by assigning positive probability everywhere non-detecting stations had
  low detection efficiency. To replicate previous behaviour, set
  `integration_method = "additive"`.

## New Features

* **Integration method parameter**: New `integration_method` parameter in
  `calculate_fish_positions()` with three options:
  - `"subtractive"` (new default): Detection field minus non-detection penalty,
    clamped to 0. Produces detection-anchored positions.
  - `"multiplicative"`: Detection field scaled down proportionally by
    non-detection evidence. Smoother penalty than subtractive.
  - `"additive"`: Original WADE formula preserved for backward compatibility.

* For subtractive and multiplicative methods, `detection_weight` is no longer
  used (detection is always the base); only `non_detection_weight` controls
  non-detection penalty strength (0-1).

* `integrated_prob` values are now rescaled to [0, 1] per fish/time period
  after integration, ensuring consistent probability ranges across all
  integration methods.

---

# positionR 1.3.0

## New Features

* **Depth-dependent behavioural state bias**: New `depth_state_bias` parameter in `simulate_fish_tracks()` biases state transition probabilities based on fish depth relative to species-specific depth thresholds, producing more ecologically realistic movement patterns.
* **Goal-directed movement (BCRW)**: New `goal_locations`, `goal_bias`, and `goal_tolerance` parameters enable biased correlated random walks where fish navigate toward specified destinations with tunable directedness.
* **`path_color` parameter**: New option in `plot_fish_tracks()` to override automatic colouring with a single colour for all paths.
* **New vignette**: `Simulation_Goal_Directed.R` demonstrates goal-directed movement across a range of bias strengths.

---

# positionR 1.2.0

## New Features

### Barrier Masking for Realistic Detection Modeling

* Added barrier masking functionality to prevent unrealistic detections through land obstacles (islands, peninsulas, shorelines)
* New `include_barriers` parameter available in four core functions:
  - `calculate_station_distances()`: Automatically detects line-of-sight barriers using ray-tracing algorithm
  - `calculate_detection_system()`: Incorporates barriers into system-wide coverage calculations
  - `simulate_fish_tracks()`: Prevents simulated detections through land obstacles
  - `calculate_fish_positions()`: Accounts for barriers in both detection and non-detection probabilities

### How Barrier Masking Works

* Line-of-sight analysis identifies when direct path between receiver and location crosses land
* Detection efficiency automatically set to 0 where barriers present
* Works with both static DE mode (pre-computed values) and temporal DE mode (on-the-fly calculation)
* Generates `crosses_barrier` column in station distances data for downstream use

## Improvements

* Enhanced realism in detection simulations by accounting for physical geography
* More accurate position estimates that respect landscape features
* Updated vignettes with barrier masking examples and visualizations:
  - Array Design & Simulation vignette: barrier field visualization and detection masking
  - WADE Positioning vignette: barrier-aware positioning workflow

## Documentation

* Added barrier masking sections to both major vignettes
* New visualizations showing barrier fields from receiver perspectives
* Documentation for `include_barriers` parameter across all affected functions
* Alert boxes explaining barrier masking benefits and requirements

---

# positionR 1.1.0

## New Features

### Interactive Array Design Application

* Added `run_array_design()` function to launch interactive Shiny app for acoustic receiver array design and optimization
* App includes four main modules:
  - **Array Design**: Interactive map for manual receiver placement and regular array generation
  - **Range Modeling**: Depth-dependent detection efficiency model fitting
  - **Array Performance**: Quantitative coverage assessment and sampling bias analysis
  - **Instructions**: Comprehensive methodology guide and best practices

### App Features

* Load custom raster files (.tif, .grd, .asc) or use included sample depth data
* Place receivers manually via click-to-add or generate regular arrays
* Calculate cost-weighted distances for realistic coverage assessment
* Visualize detection efficiency surfaces using fitted models
* Export receiver coordinates as CSV for field deployment
* Analyze depth sampling bias and coverage thresholds (>5%, >25%, >50%, >75% DE)
* Real-time performance metrics with color-coded indicators

## Improvements

* Integrated array design app with existing package functions (`generate_regular_points()`, `create_logistic_curve_depth()`)
* Enhanced namespace handling to avoid conflicts between raster and dplyr packages
* Added comprehensive documentation for array design workflow

## Dependencies

* Added Shiny app dependencies to Suggests: shiny, shinydashboard, plotly, DT

---

# positionR 1.0.0

* Initial CRAN-ready release
* Comprehensive tools for acoustic telemetry array design, fish movement simulation, and WADE positioning methodology
