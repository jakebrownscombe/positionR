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
