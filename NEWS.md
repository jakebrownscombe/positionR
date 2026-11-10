# positionR 2.0.0

- positionR 2.0 replaces the single **positionR** 1.x package with four packages:
  **detqc** reads, checks and filters detections; **detfield** builds detection fields and
  positions animals; **arraydesign** designs arrays and tests the performance of tracking
  systems; **arraydesignR** runs **arraydesign** as a point-and-click application.
- This package, **positionR** 2.0.0, holds no analysis code. Installing it installs the four
  packages, and `positionR_packages()` lists them with their installed versions and sites. It
  does not attach them.
- **detfield** positions each animal as a location surface per time bin, an hour by default.
  The surface is the probability of each water cell, by default from a binomial likelihood over
  every receiver that was listening. **positionR** 1.x gave weighted average detection
  efficiency (WADE) positions.
- **positionR** 1.7.0, the last 1.x release, stays installable from the tag `v1.7.0`:
  `remotes::install_github("jakebrownscombe/positionR@v1.7.0")`.

## Where the 1.x functions went

Every exported function of **positionR** 1.7.0 is in the table once. "No 2.0 function" means the
four packages export nothing that does that job.

| Job in 1.x | 1.x functions | In 2.0 |
|---|---|---|
| Place receivers on the water | `generate_regular_points()`, `generate_exact_regular_points()`, `generate_spaced_points()`, `generate_random_points()`, `plot_points_on_input()` | `place_receivers()` and its `plot()` method (**arraydesign**) |
| Detection range from depth | `create_logistic_curve_depth()` | `de_logistic()`; a fitted range model through `range_test_counts()` and `de_model()` (**detfield**) |
| Distances between receivers and cells | `calculate_station_distances()` | `water_distances()`, measured through water (**detfield**) |
| Detection probability of the array | `calculate_detection_system()` | `de_field()` and `system_de()` (**detfield**); `coverage_summary()` (**arraydesign**) |
| Simulate tracks and detections | `simulate_fish_tracks()`, `plot_fish_tracks()` | `simulate_tracks()`, `simulate_transmissions()`, `simulate_detections()`, `crw_model()`, `movement_model()`, `track_control()` and the `plot()` method for tracks (**arraydesign**) |
| Detection rates of simulated animals | `calculate_detection_summaries()`, `print_detection_summary()`, `plot_detection_rates()`, `analyze_detection_performance()` | Partly: the print of `simulate_detections()` gives the share of transmissions heard (**arraydesign**), and the Simulation tab of `run_app()` (**arraydesignR**) shows the rate of each animal. No 2.0 function returns or plots rates per animal |
| Score positions against known tracks | `analyze_positioning_performance()` | `score_surfaces()` and `evaluate_array()` (**arraydesign**) |
| Interactive array design | `run_array_design()` | `run_app()` (**arraydesignR**) |
| Prepare detections for positioning | `prepare_detection_data_for_wade()` | `read_detections()`, `read_deployments()`, the checks and filters, and `as_detfield_input()` (**detqc**); `bin_detections()` (**detfield**) |
| Position animals | `calculate_fish_positions()`, `plot_fish_positions()` | `location_surface()` and its `plot()` method, `credible_region()` (**detfield**) |
| Sample points from position surfaces | `sample_points_from_probabilities()` | `sample_locations()` (**detfield**) |
| Sample points from the array's detection probability | `sample_points_from_system_de()` | `sample_available()` (**detfield**): availability cells from the whole detection surface, in proportion to the probability that an hour is detected from each cell, or uniformly above a threshold, per animal and hour |
| Sample points from track grid cells | `sample_points_from_track_grid()` | No 2.0 function |
| Space use over a period | `calculate_space_use()`, `plot_space_use()`, `compare_space_use_methods()`, `plot_space_use_method_comparison()`, `compare_space_use_thresholds()`, `generate_random_points_in_space_use()`, `plot_space_use_points()`, `summarize_space_use_points()` | No 2.0 function. `credible_region()` (**detfield**) gives the region of one animal and hour |
| Presence and absence points | `generate_space_use_absences()`, `generate_space_use_absences_simple()`, `generate_simple_raster_absences()`, `test_raster_absence_generation()`, `create_presence_absence_dataset()`, `plot_presence_absence_space_use()` | Partly: presence points from position surfaces with `sample_locations()`, availability points with `sample_available()`, stacked into one table by `used_available()` (**detfield**). No 2.0 function draws absences from a space-use region or plots them |
| Habitat and depth selection | `analyze_comparative_habitat_selection()`, `print_comparative_summary()`, `plot_comparative_results()`, `plot_depth_selection()`, `plot_depth_comparison()` | Partly: `used_available()` (**detfield**) builds the table a conditional logistic regression takes and fits no model; the article "Habitat selection: presence and availability from positions" fits one and shows that position draws underestimate a preference for shallow water. `habitat_selection()` and `available_steps()` (**arraydesign**) simulate animals with a known preference and draw their control steps. No 2.0 function summarises or plots selection |
| Weight receivers by information | `calculate_station_de_characteristics()`, `apply_information_weighting()`, `print_station_characteristics()` | No 2.0 function |
| Behavioural states and temperature | `analyze_behavioral_temperature()`, `plot_behavioral_states()` | No 2.0 function |
| Example data | `stoney_fish_detections`, `stoney_rx_deploy`, walleye detections and receivers from Stoney Lake | No 2.0 dataset holds these detections. The 2.0 examples are simulated studies of Stoney Lake |
