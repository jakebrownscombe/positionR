# <img src="man/figures/logo.png" align="left" height="120" alt="" />


#positionR
Acoustic Telemetry System Design and Fish Movement Simulation
positionR is an R package for designing acoustic telemetry systems and simulating fish movement with realistic detection patterns. The package provides tools for receiver placement optimization, detection efficiency modeling, and fish track simulation with detection events.
Features

Point Generation: Create receiver station layouts with regular, spaced, or random patterns
Detection Modeling: Build realistic detection efficiency models that vary with distance and depth
System Analysis: Calculate detection probabilities across spatial areas for system optimization
Movement Simulation: Simulate realistic fish movement using correlated random walks
Detection Simulation: Model acoustic detection events based on receiver proximity and environmental factors
Visualization: Comprehensive plotting functions for tracks, detections, and system performance

Installation
Install the development version from GitHub:
r# Install devtools if you haven't already
install.packages("devtools")

## Install positionR
devtools::install_github("yourusername/positionR", build_vignettes = TRUE)
Quick Start
rlibrary(positionR)

## Load example depth data (included with package)
data(depth_raster)

## 1. Generate receiver station locations
stations <- generate_random_points(depth_raster, n_points = 6, seed = 123)

## 2. Calculate distances between stations and all raster cells
distances <- calculate_station_distances(depth_raster, stations, max_distance = 500)

## 3. Create detection efficiency model
de_model <- create_logistic_curve_depth(
  min_depth = 2, max_depth = 30,
  d50_min_depth = 50, d95_min_depth = 20,
  d50_max_depth = 150, d95_max_depth = 60,
  plot = FALSE
)

## 4. Calculate system detection probabilities
detection_system <- calculate_detection_system(
  distance_frame = distances,
  receiver_frame = stations,
  model = de_model$log_model,
  output_type = "cumulative"
)

## 5. Simulate fish tracks with detections
fish_tracks <- simulate_fish_tracks(
  raster = depth_raster,
  detection_system = detection_system$data,
  station_distances = distances,
  n_paths = 5,
  n_steps = 100
)

## 6. Visualize results
plot_fish_tracks(fish_tracks, depth_raster, stations)

## 7. Analyze detection performance
analyze_detection_performance(fish_tracks)



##Core Functions


##Receiver Placement

generate_exact_regular_points() - Generate specific number of regularly spaced points
generate_spaced_points() - Generate points with defined spacing
generate_random_points() - Generate randomly distributed points

##Detection Modeling

create_logistic_curve_depth() - Build detection efficiency models varying with distance and depth
calculate_station_distances() - Calculate cost-weighted distances from receivers to all locations
calculate_detection_system() - Compute system-wide detection probabilities

##Movement Simulation

simulate_fish_tracks() - Generate realistic fish movement paths with detection events
plot_fish_tracks() - Visualize tracks, detections, and receiver performance

##Analysis

calculate_detection_summaries() - Compute detection rate statistics
analyze_detection_performance() - Comprehensive detection analysis with plots

##Workflow
The typical workflow follows these steps:

Design receiver array using point generation functions
Model detection efficiency with distance/depth relationships
Calculate system coverage to evaluate detection probabilities
Simulate fish movement with realistic detection events
Analyze performance to optimize system design

Example Output
Detection System Coverage
Show Image
Fish Tracks with Detections
Show Image
Detection Performance Analysis
Show Image
Key Applications

System Design: Optimize receiver placement for maximum detection coverage
Performance Evaluation: Assess detection efficiency of existing arrays
Movement Ecology: Study fish movement patterns and habitat use
Detection Modeling: Understand factors affecting acoustic detection range
Simulation Studies: Test system performance under different scenarios

##Vignette
For detailed examples and tutorials, see the package vignette:
rvignette("positionR_tutorial", package = "positionR")
Dependencies
The package builds on several key R packages:

raster - Spatial raster data handling
sf - Simple features for spatial vector data
ggplot2 - Data visualization
dplyr - Data manipulation
gdistance - Cost-distance calculations
circular - Circular statistics for movement modeling

##Citation
If you use positionR in your research, please cite:
Brownscombe, J.W. (2025). positionR: Acoustic Telemetry System Design and Fish Movement Simulation. 
R package version [version]. https://github.com/jakebrownscombe/positionR

##Contributing
Contributions are welcome! Please feel free to submit a Pull Request.


Contact

Author: Dr. Jacob Brownscombe
Email: jakebrownscombe@gmail.com
GitHub: @jakebrownscombe

