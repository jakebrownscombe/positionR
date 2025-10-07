library(shiny)
library(plotly)
library(raster)
library(DT)
library(ggplot2)
library(sf)
library(methods)
library(gdistance)
library(dplyr)
library(tidyr)
library(positionR)  # Load package functions

# Helper function - local implementation of calculate_station_distances
# with dplyr:: namespace fixes for compatibility
calculate_station_distances_app <- function(raster, receiver_frame, max_distance = NULL, station_col = "point_id") {

  # Validate inputs
  if (!station_col %in% names(receiver_frame)) {
    stop("Column '", station_col, "' not found in receiver_frame.")
  }

  # Get station IDs from specified column
  station_ids <- receiver_frame[[station_col]]
  if (any(duplicated(station_ids))) {
    stop("Station IDs must be unique.")
  }

  # Create uniform cost surface
  cost_raster <- raster
  cost_raster[!is.na(cost_raster)] <- 1

  # Create transition matrix
  tr <- transition(cost_raster, transitionFunction = mean, directions = 8)
  tr <- geoCorrection(tr, type = "c")

  # Convert points to SpatialPoints
  points_sp <- as(receiver_frame, "Spatial")

  # Get raster cell coordinates
  raster_coords <- coordinates(cost_raster)
  cell_numbers <- 1:ncell(cost_raster)
  valid_cells <- which(!is.na(values(cost_raster)))

  # Initialize results dataframe
  results_df <- data.frame(
    cell_id = cell_numbers[valid_cells],
    x = raster_coords[valid_cells, 1],
    y = raster_coords[valid_cells, 2],
    raster_value = values(raster)[valid_cells]
  )

  # Calculate cost distances from each station
  for (i in 1:nrow(receiver_frame)) {
    current_station_id <- station_ids[i]

    # Get single point
    single_point <- points_sp[i, ]
    station_coords <- coordinates(single_point)

    # Calculate cost distance
    cost_dist_raster <- accCost(tr, single_point)
    cost_distances <- values(cost_dist_raster)[valid_cells]

    # Calculate straight-line distances
    cell_coords <- raster_coords[valid_cells, ]
    straight_distances <- sqrt((cell_coords[,1] - station_coords[1])^2 +
                                 (cell_coords[,2] - station_coords[2])^2)

    # Replace infinite values with NA
    cost_distances[is.infinite(cost_distances)] <- NA

    # Apply maximum distance filter
    if (!is.null(max_distance)) {
      cost_distances[cost_distances > max_distance] <- NA
      straight_distances[straight_distances > max_distance] <- NA
    }

    # Add to results
    safe_station_id <- gsub("[^[:alnum:]_]", "_", as.character(current_station_id))
    results_df[[paste0("cost_dist_station_", safe_station_id)]] <- cost_distances
    results_df[[paste0("straight_dist_station_", safe_station_id)]] <- straight_distances
  }

  # Convert to long format
  cost_cols <- grep("^cost_dist_station_", names(results_df), value = TRUE)
  straight_cols <- grep("^straight_dist_station_", names(results_df), value = TRUE)

  # Create mapping
  station_mapping <- data.frame(
    safe_name = gsub("[^[:alnum:]_]", "_", as.character(station_ids)),
    original_id = station_ids,
    stringsAsFactors = FALSE
  )

  # Pivot to long format - use dplyr:: to avoid namespace conflicts with raster package
  cost_long <- results_df %>%
    dplyr::select(cell_id, x, y, raster_value, all_of(cost_cols)) %>%
    pivot_longer(cols = all_of(cost_cols), names_to = "station_col", values_to = "cost_distance") %>%
    mutate(safe_name = gsub("cost_dist_station_", "", station_col)) %>%
    left_join(station_mapping, by = "safe_name") %>%
    rename(station_no = original_id) %>%
    dplyr::select(-station_col, -safe_name)

  straight_long <- results_df %>%
    dplyr::select(cell_id, all_of(straight_cols)) %>%
    pivot_longer(cols = all_of(straight_cols), names_to = "station_col", values_to = "straight_distance") %>%
    mutate(safe_name = gsub("straight_dist_station_", "", station_col)) %>%
    left_join(station_mapping, by = "safe_name") %>%
    rename(station_no = original_id) %>%
    dplyr::select(-station_col, -safe_name)

  # Combine and calculate tortuosity
  final_df <- cost_long %>%
    left_join(straight_long, by = c("cell_id", "station_no")) %>%
    mutate(tortuosity = cost_distance / straight_distance) %>%
    filter(!is.na(cost_distance)) %>%
    dplyr::select(cell_id, x, y, raster_value, station_no, cost_distance, straight_distance, tortuosity) %>%
    arrange(station_no, cell_id)

  return(final_df)
}

server <- function(input, output, session) {

  # Reactive values to store data
  values <- reactiveValues(
    raster_data = NULL,
    raster_df = NULL,
    receivers = data.frame(
      id = integer(),
      x = numeric(),
      y = numeric(),
      depth = numeric(),
      type = character(),  # "manual" or "regular"
      stringsAsFactors = FALSE
    ),
    next_id = 1,
    distances_calculated = FALSE,
    station_distances = NULL,
    detection_model = NULL
  )

  # Auto-fit range model on app startup with default settings
  observe({
    # Only run once when session starts
    if (is.null(values$detection_model)) {
      tryCatch({
        values$detection_model <- positionR::create_logistic_curve_depth(
          min_depth = 1,        # Default from UI
          max_depth = 35,       # Default from UI
          d50_min_depth = 400,  # Default from UI
          d95_min_depth = 800,  # Default from UI
          d50_max_depth = 750,  # Default from UI
          d95_max_depth = 1500, # Default from UI
          plot = FALSE
        )

        showNotification("Default detection efficiency model loaded!", type = "message", duration = 3)

      }, error = function(e) {
        showNotification(paste("Error loading default model:", e$message), type = "warning", duration = 5)
      })
    }
  })

  # Load raster function with CRS handling
  load_raster <- function(file_path, file_name = "raster") {
    tryCatch({
      # Load the raster
      r <- raster(file_path)
      original_crs <- crs(r)

      print(paste("Raster loaded - dimensions:", nrow(r), "x", ncol(r)))
      print(paste("Original CRS:", as.character(original_crs)))

      # Check if CRS is geographic (lat/lon)
      is_geographic <- !is.na(original_crs) && grepl("\\+proj=longlat", as.character(original_crs))

      # If geographic, project to appropriate UTM zone for distance calculations
      if (is_geographic) {
        # Calculate center point to determine UTM zone
        center_lon <- mean(c(extent(r)@xmin, extent(r)@xmax))
        center_lat <- mean(c(extent(r)@ymin, extent(r)@ymax))

        # Calculate UTM zone
        utm_zone <- floor((center_lon + 180) / 6) + 1

        # Determine hemisphere
        hemisphere <- if (center_lat >= 0) "north" else "south"

        # Create UTM CRS string
        utm_crs <- paste0("+proj=utm +zone=", utm_zone,
                          if (hemisphere == "south") " +south" else "",
                          " +datum=WGS84 +units=m +no_defs")

        print(paste("Converting from Geographic to UTM Zone", utm_zone, hemisphere))
        showNotification(paste("Converting coordinates from Lat/Lon to UTM Zone", utm_zone, hemisphere, "for distance calculations..."),
                         type = "message", duration = 5)

        # Project the raster
        r_proj <- projectRaster(r, crs = utm_crs, method = "bilinear")

        # Store both original and projected
        values$raster_data <- r_proj  # Use projected for calculations
        values$original_raster <- r   # Keep original for reference
        values$original_crs <- original_crs
        values$projected_crs <- utm_crs

        r <- r_proj  # Use projected raster for data frame creation

      } else {
        # Already projected, use as-is
        values$raster_data <- r
        values$original_raster <- r
        values$original_crs <- original_crs
        values$projected_crs <- original_crs

        print("Raster already in projected coordinates")
      }

      # Convert to data frame step by step
      coords <- coordinates(r)
      vals <- values(r)

      # Create data frame manually
      raster_df <- data.frame(
        x = coords[, 1],
        y = coords[, 2],
        value = as.numeric(vals)
      )

      print(paste("DataFrame created - dimensions:", nrow(raster_df), "x", ncol(raster_df)))

      # Remove NA values
      raster_df <- raster_df[!is.na(raster_df$value), ]

      if (nrow(raster_df) == 0) {
        stop("No valid raster data found")
      }

      values$raster_df <- raster_df

      # Show appropriate notification
      crs_info <- if (is_geographic) paste("(converted to UTM Zone", utm_zone, hemisphere, ")") else "(projected coordinates)"
      showNotification(paste(file_name, "loaded! Size:", nrow(raster_df), "pixels", crs_info))

    }, error = function(e) {
      showNotification(paste("Error loading", file_name, ":", as.character(e)), type = "error")
      print(paste("Full error:", e))
    })
  }

  # Load default raster on startup - check for files in data/ subdirectory
  observe({
    # Only load once when app starts
    if (is.null(values$raster_df)) {
      # Try data/ subdirectory files first
      app_dir <- system.file("shiny-apps/array_design", package = "positionR")
      depth_raster_path <- file.path(app_dir, "data", "depth_raster.tif")
      depth_path <- file.path(app_dir, "data", "depth.tif")

      if (file.exists(depth_raster_path)) {
        load_raster(depth_raster_path, "Default depth raster")
      } else if (file.exists(depth_path)) {
        load_raster(depth_path, "Default depth raster")
      }
    }
  })

  # Load raster file from user upload
  observeEvent(input$raster_file, {
    req(input$raster_file)
    load_raster(input$raster_file$datapath, "User raster")
  })

  # Generate regular points - uses positionR package function
  observeEvent(input$generate_regular, {
    req(values$raster_data)

    tryCatch({
      # Generate regular points using positionR function
      regular_points_sf <- positionR::generate_regular_points(
        values$raster_data,
        n_points = input$n_regular_points,
        seed = input$array_seed
      )

      # Convert to data frame and add to receivers
      n_generated <- nrow(regular_points_sf)
      if (n_generated > 0) {
        new_receivers <- data.frame(
          id = values$next_id:(values$next_id + n_generated - 1),
          x = round(regular_points_sf$x, 2),
          y = round(regular_points_sf$y, 2),
          depth = round(regular_points_sf$raster_value, 2),
          type = "regular",
          stringsAsFactors = FALSE
        )

        values$receivers <- rbind(values$receivers, new_receivers)
        values$next_id <- values$next_id + n_generated

        showNotification(paste(n_generated, "regular receivers generated"))
      } else {
        showNotification("No points could be generated in the raster area", type = "warning")
      }

    }, error = function(e) {
      showNotification(paste("Error generating regular points:", e$message), type = "error")
      print(paste("Regular point error:", e))
    })
  })

  # Calculate distances - uses local helper function with namespace fixes
  observeEvent(input$calculate_distances, {
    req(values$raster_data)
    req(nrow(values$receivers) > 0)

    tryCatch({
      showNotification("Calculating distances... This may take several minutes.", type = "message")

      # Convert receivers to sf object for distance calculation
      receivers_sf <- st_as_sf(values$receivers, coords = c("x", "y"), crs = st_crs(values$raster_data))
      receivers_sf$point_id <- values$receivers$id

      # Calculate distances using local helper function
      distances <- calculate_station_distances_app(
        raster = values$raster_data,
        receiver_frame = receivers_sf,
        max_distance = input$max_distance,
        station_col = "point_id"
      )

      values$station_distances <- distances
      values$distances_calculated <- TRUE

      # Automatically turn on show distance fields checkbox
      updateCheckboxInput(session, "show_distances", value = TRUE)

      showNotification(paste("Distance calculations complete!", nrow(distances), "station-cell combinations calculated"))

    }, error = function(e) {
      showNotification(paste("Error calculating distances:", e$message), type = "error")
      print(paste("Distance calculation error:", e))
    })
  })

  # Create the interactive plot
  output$raster_plot <- renderPlotly({
    req(values$raster_df)

    # Create basic plotly plot
    p <- plot_ly(source = "raster_plot")

    # Show detection efficiency if toggle is enabled and both model and distances are available
    if (input$show_detection_efficiency && !is.null(values$detection_model) && !is.null(values$station_distances)) {
      # Calculate detection efficiency for each pixel using the fitted model
      de_df <- values$station_distances %>%
        filter(!is.na(cost_distance)) %>%
        mutate(de_individual = mapply(values$detection_model$func, cost_distance, raster_value)) %>%
        group_by(x, y, raster_value) %>%
        # Calculate system-wide detection efficiency (1 - product of all non-detections)
        summarise(detection_efficiency = 1 - prod(1 - de_individual, na.rm = TRUE), .groups = 'drop')

      # Magma colorscale for detection efficiency (purple = low, yellow = high)
      de_colorscale <- list(
        c(0, "#000004"),   # Dark purple (low detection efficiency)
        c(0.2, "#3B0F70"),
        c(0.4, "#8E1538"),
        c(0.6, "#FD8861"),
        c(0.8, "#FEC488"),
        c(1, "#FCFFA4")    # Light yellow (high detection efficiency)
      )

      p <- p %>%
        add_heatmap(
          data = de_df,
          x = ~x,
          y = ~y,
          z = ~detection_efficiency,
          colorscale = de_colorscale,
          hovertemplate = "X: %{x}<br>Y: %{y}<br>Depth: %{customdata}<br>Detection Efficiency: %{z:.2f}<extra></extra>",
          customdata = ~raster_value,
          name = "Detection Efficiency"
        )
    } else if (input$show_distances && !is.null(values$station_distances)) {
      # Create distance field visualization - use minimum distance across all receivers
      distance_df <- values$station_distances %>%
        filter(!is.na(cost_distance)) %>%
        group_by(x, y) %>%
        summarise(min_distance = min(cost_distance, na.rm = TRUE), .groups = 'drop')

      # Red-to-yellow distance colorscale (red = near, yellow = far)
      distance_colorscale <- list(
        c(0, "red"),
        c(0.3, "orange"),
        c(0.6, "yellow"),
        c(1, "lightyellow")
      )

      p <- p %>%
        add_heatmap(
          data = distance_df,
          x = ~x,
          y = ~y,
          z = ~min_distance,
          colorscale = distance_colorscale,
          hovertemplate = "X: %{x}<br>Y: %{y}<br>Min Distance: %{z:.0f}m<extra></extra>",
          name = "Distance Field"
        )
    } else {
      # Show depth data as before
      blue_colorscale <- list(
        c(0, "darkblue"),     # deepest water (most negative/lowest values)
        c(0.2, "blue"),
        c(0.4, "mediumblue"),
        c(0.6, "dodgerblue"),
        c(0.8, "deepskyblue"),
        c(1, "skyblue")       # shallowest water (higher values)
      )

      p <- p %>%
        add_heatmap(
          data = values$raster_df,
          x = ~x,
          y = ~y,
          z = ~value,
          colorscale = blue_colorscale,
          hovertemplate = "X: %{x}<br>Y: %{y}<br>Depth: %{z}<extra></extra>",
          name = "Depth"
        )
    }

    # Add receiver points if any exist - different colors for manual vs regular
    if (nrow(values$receivers) > 0) {
      # Manual receivers (red)
      manual_receivers <- values$receivers[values$receivers$type == "manual" | is.na(values$receivers$type), ]
      if (nrow(manual_receivers) > 0) {
        p <- p %>%
          add_markers(
            data = manual_receivers,
            x = ~x,
            y = ~y,
            marker = list(
              size = 10,
              color = "red",
              line = list(color = "white", width = 1.5)
            ),
            hovertemplate = "Manual Receiver<br>ID: %{text}<br>X: %{x}<br>Y: %{y}<br>Depth: %{customdata}<extra></extra>",
            text = ~id,
            customdata = ~depth,
            name = "Manual"
          )
      }

      # Regular receivers (orange)
      regular_receivers <- values$receivers[values$receivers$type == "regular", ]
      if (nrow(regular_receivers) > 0) {
        p <- p %>%
          add_markers(
            data = regular_receivers,
            x = ~x,
            y = ~y,
            marker = list(
              size = 10,
              color = "orange",
              line = list(color = "white", width = 1.5)
            ),
            hovertemplate = "Regular Array<br>ID: %{text}<br>X: %{x}<br>Y: %{y}<br>Depth: %{customdata}<extra></extra>",
            text = ~id,
            customdata = ~depth,
            name = "Regular"
          )
      }
    }

    # Calculate aspect ratio to prevent stretching
    x_range <- range(values$raster_df$x)
    y_range <- range(values$raster_df$y)
    x_span <- diff(x_range)
    y_span <- diff(y_range)

    # Dynamic title based on what's being shown
    plot_title <- if (input$show_detection_efficiency && !is.null(values$detection_model) && !is.null(values$station_distances)) {
      "Detection Efficiency - Click to Add/Remove Receivers"
    } else if (input$show_distances && !is.null(values$station_distances)) {
      "Distance Fields - Click to Add/Remove Receivers"
    } else {
      "Depth Map - Click to Add/Remove Receivers"
    }

    # Simple layout with proper aspect ratio
    p %>%
      layout(
        title = plot_title,
        xaxis = list(
          title = "X Coordinate",
          scaleanchor = "y",
          scaleratio = 1
        ),
        yaxis = list(
          title = "Y Coordinate",
          scaleanchor = "x"
        )
      ) %>%
      config(displayModeBar = TRUE, displaylogo = FALSE)
  })

  # Handle click events to add/remove receivers
  observeEvent(event_data("plotly_click", source = "raster_plot"), {
    click_data <- event_data("plotly_click", source = "raster_plot")

    if (!is.null(click_data) && !is.null(click_data$x) && !is.null(click_data$y)) {

      # Check if click was on a receiver (curveNumber > 0 means markers layer)
      if (!is.null(click_data$curveNumber) && click_data$curveNumber > 0) {
        # Debug: print click info
        print(paste("Clicked curveNumber:", click_data$curveNumber, "pointNumber:", click_data$pointNumber))

        # Clicked on a receiver - determine which type and remove it
        point_index <- click_data$pointNumber + 1

        # Get the actual receiver based on which layer was clicked
        manual_receivers <- values$receivers[values$receivers$type == "manual" | is.na(values$receivers$type), ]
        regular_receivers <- values$receivers[values$receivers$type == "regular", ]

        # Determine layer mapping based on what exists
        manual_layer <- if (nrow(manual_receivers) > 0) 1 else NULL
        regular_layer <- if (nrow(regular_receivers) > 0) {
          if (is.null(manual_layer)) 1 else 2
        } else NULL

        if (click_data$curveNumber == manual_layer && !is.null(manual_layer)) {
          # Manual receivers
          if (point_index <= nrow(manual_receivers)) {
            receiver_to_remove <- manual_receivers[point_index, ]
            receiver_row <- which(values$receivers$id == receiver_to_remove$id)
            if (length(receiver_row) > 0) {
              removed_id <- values$receivers$id[receiver_row[1]]
              values$receivers <- values$receivers[-receiver_row[1], ]
              showNotification(paste("Manual receiver", removed_id, "removed"))
            }
          }
        } else if (click_data$curveNumber == regular_layer && !is.null(regular_layer)) {
          # Regular receivers
          if (point_index <= nrow(regular_receivers)) {
            receiver_to_remove <- regular_receivers[point_index, ]
            receiver_row <- which(values$receivers$id == receiver_to_remove$id)
            if (length(receiver_row) > 0) {
              removed_id <- values$receivers$id[receiver_row[1]]
              values$receivers <- values$receivers[-receiver_row[1], ]
              showNotification(paste("Regular receiver", removed_id, "removed"))
            }
          }
        }
      } else {
        # Clicked on raster - add new receiver with depth extraction
        click_x <- round(as.numeric(click_data$x), 2)
        click_y <- round(as.numeric(click_data$y), 2)

        # Extract depth value from raster at click location
        click_depth <- NA
        if (!is.null(values$raster_data)) {
          tryCatch({
            # Use extract function to get raster value at clicked coordinates
            click_point <- cbind(click_x, click_y)
            click_depth <- raster::extract(values$raster_data, click_point)
            click_depth <- round(as.numeric(click_depth), 2)
          }, error = function(e) {
            print(paste("Error extracting depth:", e$message))
          })
        }

        new_receiver <- data.frame(
          id = values$next_id,
          x = click_x,
          y = click_y,
          depth = click_depth,
          type = "manual",
          stringsAsFactors = FALSE
        )

        values$receivers <- rbind(values$receivers, new_receiver)
        values$next_id <- values$next_id + 1

        showNotification(paste("Receiver", new_receiver$id, "added at", new_receiver$x, ",", new_receiver$y, "- Depth:", new_receiver$depth))
      }
    }
  })

  # Note: Drag functionality would require additional implementation
  # with plotly shapes or custom event handling

  # Clear all receivers - show confirmation dialog
  observeEvent(input$clear_receivers, {
    showModal(modalDialog(
      title = "Clear All Receivers",
      p("Are you sure you want to clear all receivers?"),
      p("This action cannot be undone and will also reset all distance calculations."),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_clear", "Yes, Clear All", class = "btn-danger")
      )
    ))
  })

  # Handle confirmed clear action
  observeEvent(input$confirm_clear, {
    values$receivers <- data.frame(
      id = integer(),
      x = numeric(),
      y = numeric(),
      depth = numeric(),
      type = character(),
      stringsAsFactors = FALSE
    )
    values$next_id <- 1

    # Reset distance-related data
    values$distances_calculated <- FALSE
    values$station_distances <- NULL

    # Automatically turn off both visualization checkboxes
    updateCheckboxInput(session, "show_distances", value = FALSE)
    updateCheckboxInput(session, "show_detection_efficiency", value = FALSE)

    # Close the modal and show notification
    removeModal()
    showNotification("All receivers cleared", type = "message")
  })

  # Display raster information
  output$raster_info <- renderText({
    if (is.null(values$raster_data)) {
      "No raster loaded"
    } else {
      # Build info text
      info_lines <- c(
        "Raster Information:",
        paste("Dimensions:", nrow(values$raster_data), "x", ncol(values$raster_data)),
        paste("Resolution:", round(res(values$raster_data)[1], 2), "x", round(res(values$raster_data)[2], 2), "m"),
        paste("Extent:", paste(round(extent(values$raster_data)[], 0), collapse = ", "))
      )

      # Add CRS information
      if (!is.null(values$original_crs) && !is.null(values$projected_crs)) {
        if (as.character(values$original_crs) != as.character(values$projected_crs)) {
          # Coordinate transformation occurred
          info_lines <- c(info_lines,
                          "Coordinate System: Geographic → Projected",
                          paste("Working CRS:",
                                if (grepl("utm", values$projected_crs)) {
                                  # Extract UTM zone info
                                  zone_match <- regexpr("zone=([0-9]+)", values$projected_crs)
                                  if (zone_match > 0) {
                                    zone <- regmatches(values$projected_crs, zone_match)
                                    zone <- gsub("zone=", "UTM Zone ", zone)
                                    if (grepl("south", values$projected_crs)) zone <- paste(zone, "South")
                                    zone
                                  } else "UTM"
                                } else "Projected"
                          )
          )
        } else {
          # No transformation
          crs_str <- as.character(values$original_crs)
          if (grepl("utm", crs_str)) {
            info_lines <- c(info_lines, "Coordinate System: UTM (Projected)")
          } else if (grepl("longlat", crs_str)) {
            info_lines <- c(info_lines, "Coordinate System: Geographic (Lat/Lon)")
          } else {
            info_lines <- c(info_lines, "Coordinate System: Projected")
          }
        }
      }

      paste(info_lines, collapse = "\n")
    }
  })

  # Display receiver summary
  output$receiver_summary <- renderTable({
    if (nrow(values$receivers) == 0) {
      data.frame("Summary" = "No receivers placed")
    } else {
      depth_summary <- if (all(is.na(values$receivers$depth))) {
        "No depth data"
      } else {
        paste(round(range(values$receivers$depth, na.rm = TRUE), 2), collapse = " - ")
      }

      # Count manual vs regular receivers
      manual_count <- sum(values$receivers$type == "manual" | is.na(values$receivers$type))
      regular_count <- sum(values$receivers$type == "regular", na.rm = TRUE)

      data.frame(
        "Metric" = c("Total Receivers", "Manual (Red)", "Regular (Orange)", "X Range", "Y Range", "Depth Range"),
        "Value" = c(
          nrow(values$receivers),
          manual_count,
          regular_count,
          paste(round(range(values$receivers$x), 2), collapse = " - "),
          paste(round(range(values$receivers$y), 2), collapse = " - "),
          depth_summary
        )
      )
    }
  })

  # Receiver details table (on main page)
  output$receiver_details <- DT::renderDataTable({
    if (nrow(values$receivers) == 0) {
      data.frame(
        "ID" = integer(0),
        "X Coordinate" = numeric(0),
        "Y Coordinate" = numeric(0),
        "Depth" = numeric(0),
        "Type" = character(0)
      )
    } else {
      display_data <- values$receivers
      names(display_data) <- c("ID", "X Coordinate", "Y Coordinate", "Depth", "Type")
      display_data
    }
  }, options = list(
    pageLength = 10,
    scrollX = TRUE,
    dom = 'ftp',  # Remove length menu and info
    columnDefs = list(
      list(targets = c(1, 2, 3), className = 'dt-center')  # Center align numeric columns
    )
  ))

  # Receiver data table (on data tab)
  output$receiver_table <- DT::renderDataTable({
    if (nrow(values$receivers) == 0) {
      data.frame("No receivers placed" = character(0))
    } else {
      values$receivers
    }
  }, options = list(pageLength = 15, scrollX = TRUE))

  # Download handler for receiver coordinates (triggered by action button)
  observeEvent(input$download_receivers, {
    if (nrow(values$receivers) > 0) {
      filename <- paste("receiver_coordinates_", Sys.Date(), ".csv", sep = "")
      write.csv(values$receivers, filename, row.names = FALSE)
      showNotification(paste("File saved as:", filename), type = "message")
    } else {
      showNotification("No receivers to download", type = "warning")
    }
  })

  # Range Modeling Tab - Fit Detection Model
  observeEvent(input$fit_range_model, {
    tryCatch({
      # Create logistic detection efficiency model using positionR function
      values$detection_model <- positionR::create_logistic_curve_depth(
        min_depth = input$min_depth,
        max_depth = input$max_depth,
        d50_min_depth = input$d50_min_depth,
        d95_min_depth = input$d95_min_depth,
        d50_max_depth = input$d50_max_depth,
        d95_max_depth = input$d95_max_depth,
        plot = FALSE  # Don't show plots in console for Shiny
      )

      showNotification("Detection efficiency model fitted successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error fitting model:", e$message), type = "error")
    })
  })

  # Range model plot
  output$range_model_plot <- renderPlot({
    if (is.null(values$detection_model)) {
      # Show placeholder plot
      plot(1, type = "n", xlab = "", ylab = "",
           xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
      text(0.5, 0.5, "Click 'Fit Detection Model'\nto generate visualization",
           cex = 1.2, col = "gray50")
      return()
    }

    # Create range model visualization
    library(ggplot2)

    # Generate data for different depths
    depths <- seq(input$min_depth, input$max_depth, length.out = 5)
    dist_range <- seq(0, max(input$d95_max_depth, input$d95_min_depth), by = 10)

    plot_data <- expand.grid(distance = dist_range, depth = depths)
    plot_data$detection_efficiency <- mapply(values$detection_model$func,
                                             plot_data$distance, plot_data$depth)
    plot_data$depth_label <- paste("Depth:", round(plot_data$depth, 1), "m")

    ggplot(plot_data, aes(x = distance, y = detection_efficiency, color = factor(depth))) +
      geom_line(size = 1.2) +
      scale_color_viridis_d(name = "Depth (m)", direction = -1) +
      labs(
        title = "Detection Efficiency vs Distance by Depth",
        x = "Distance (m)",
        y = "Detection Efficiency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11)
      )
  })

  # Detection efficiency prediction plot
  output$detection_efficiency_plot <- renderPlotly({
    if (is.null(values$detection_model)) {
      # Show placeholder
      p <- plot_ly() %>%
        add_annotations(
          text = "Fit detection model to view efficiency predictions",
          x = 0.5, y = 0.5, xref = "paper", yref = "paper",
          showarrow = FALSE, font = list(size = 16, color = "gray")
        ) %>%
        layout(
          title = "Detection Efficiency Surface",
          xaxis = list(title = "Distance (m)"),
          yaxis = list(title = "Depth (m)")
        )
      return(p)
    }

    # Create detection efficiency surface
    library(plotly)

    # Generate prediction grid
    dist_seq <- seq(0, max(input$d95_max_depth, input$d95_min_depth), by = 20)
    depth_seq <- seq(input$min_depth, input$max_depth, by = 1)

    grid_data <- expand.grid(distance = dist_seq, depth = depth_seq)
    grid_data$efficiency <- mapply(values$detection_model$func,
                                   grid_data$distance, grid_data$depth)

    # Create matrix for heatmap
    efficiency_matrix <- matrix(grid_data$efficiency,
                                nrow = length(dist_seq),
                                ncol = length(depth_seq))

    plot_ly(
      x = dist_seq,
      y = depth_seq,
      z = t(efficiency_matrix),
      type = "heatmap",
      colorscale = list(c(0, "#000004FF"), c(0.1, "#1B0C41FF"), c(0.2, "#4B0C6BFF"), c(0.3, "#781C6DFF"), c(0.4, "#A52C60FF"), c(0.5, "#CF4446FF"), c(0.6, "#ED6925FF"), c(0.7, "#FB9A06FF"), c(0.8, "#F7D03CFF"), c(0.9, "#FCFFA4FF"), c(1, "#FCFDBFFF")),
      alpha = 1,
      hovertemplate = "Distance: %{x}m<br>Depth: %{y}m<br>Efficiency: %{z:.2f}<extra></extra>"
    ) %>%
      layout(
        title = "Detection Efficiency Surface",
        xaxis = list(title = "Distance (m)"),
        yaxis = list(title = "Depth (m)", autorange = "reversed"),
        font = list(size = 12)
      )
  })

  # Array Performance Tab Functions and Outputs

  # Helper function to calculate detection efficiency data when needed
  get_detection_efficiency_data <- reactive({
    if (!is.null(values$detection_model) && !is.null(values$station_distances)) {
      # Calculate detection efficiency for the system
      de_df <- values$station_distances %>%
        filter(!is.na(cost_distance)) %>%
        mutate(de_individual = mapply(values$detection_model$func, cost_distance, raster_value)) %>%
        group_by(x, y, raster_value) %>%
        summarise(detection_efficiency = 1 - prod(1 - de_individual, na.rm = TRUE), .groups = 'drop')

      return(de_df)
    }
    return(NULL)
  })

  # Value boxes for key metrics
  output$mean_de <- renderValueBox({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data)) {
      mean_de <- mean(de_data$detection_efficiency, na.rm = TRUE)
      valueBox(
        value = paste0(round(mean_de * 100, 1), "%"),
        subtitle = "Mean Detection Efficiency",
        icon = icon("bullseye"),
        color = if(mean_de > 0.5) "green" else if(mean_de > 0.25) "yellow" else "red"
      )
    } else {
      valueBox(
        value = "N/A",
        subtitle = "Mean Detection Efficiency",
        icon = icon("bullseye"),
        color = "light-blue"
      )
    }
  })

  output$coverage_area <- renderValueBox({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data)) {
      total_pixels <- nrow(de_data)
      covered_pixels <- sum(de_data$detection_efficiency > 0.05, na.rm = TRUE)
      coverage_pct <- (covered_pixels / total_pixels) * 100

      valueBox(
        value = paste0(round(coverage_pct, 1), "%"),
        subtitle = "Area with >5% DE",
        icon = icon("map"),
        color = if(coverage_pct > 75) "green" else if(coverage_pct > 50) "yellow" else "red"
      )
    } else {
      valueBox(
        value = "N/A",
        subtitle = "Area with >5% DE",
        icon = icon("map"),
        color = "light-blue"
      )
    }
  })

  output$high_de_area <- renderValueBox({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data)) {
      total_pixels <- nrow(de_data)
      high_de_pixels <- sum(de_data$detection_efficiency > 0.5, na.rm = TRUE)
      high_de_pct <- (high_de_pixels / total_pixels) * 100

      valueBox(
        value = paste0(round(high_de_pct, 1), "%"),
        subtitle = "Area with >50% DE",
        icon = icon("star"),
        color = if(high_de_pct > 50) "green" else if(high_de_pct > 25) "yellow" else "red"
      )
    } else {
      valueBox(
        value = "N/A",
        subtitle = "Area with >50% DE",
        icon = icon("star"),
        color = "light-blue"
      )
    }
  })

  # Detection efficiency histogram
  output$de_histogram <- renderPlot({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data)) {
      library(ggplot2)

      # Subsample for performance if dataset is large
      if (nrow(de_data) > 5000) {
        de_sample <- de_data[sample(nrow(de_data), 5000), ]
      } else {
        de_sample <- de_data
      }

      # Create dataframe for threshold lines with labels
      threshold_data <- data.frame(
        x = c(0.05, 0.25, 0.5, 0.75),
        label = c("5% DE", "25% DE", "50% DE", "75% DE"),
        color = c("red", "darkorange", "gold", "darkgreen")
      )

      ggplot(de_sample, aes(x = detection_efficiency)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
        geom_vline(data = threshold_data, aes(xintercept = x, color = label),
                   linetype = "dashed", size = 1) +
        scale_color_manual(
          name = "Coverage Thresholds",
          values = c("5% DE" = "red", "25% DE" = "darkorange", "50% DE" = "gold", "75% DE" = "darkgreen")
        ) +
        labs(
          title = "Distribution of Detection Efficiency",
          x = "Detection Efficiency",
          y = "Frequency"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 10),
          legend.position = "right",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8)
        )
    } else {
      # Placeholder plot
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Calculate distances and fit model\nto view performance metrics",
                 size = 4, color = "gray50") +
        theme_void()
    }
  })

  # Coverage thresholds table
  output$coverage_thresholds <- renderTable({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data)) {
      total_pixels <- nrow(de_data)

      thresholds <- c(0.05, 0.25, 0.5, 0.75)
      results <- data.frame(
        Threshold = paste0(">", thresholds * 100, "% DE"),
        Pixels = sapply(thresholds, function(t) sum(de_data$detection_efficiency > t, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      results$`Area (%)` <- round((results$Pixels / total_pixels) * 100, 1)

      results[, c("Threshold", "Area (%)", "Pixels")]
    } else {
      data.frame(
        Threshold = c(">5% DE", ">25% DE", ">50% DE", ">75% DE"),
        `Area (%)` = c("--", "--", "--", "--"),
        Pixels = c("--", "--", "--", "--"),
        stringsAsFactors = FALSE
      )
    }
  }, striped = TRUE, hover = TRUE)

  # Depth coverage comparison
  output$depth_comparison <- renderPlot({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data) && !is.null(values$raster_df)) {
      library(ggplot2)
      library(patchwork)

      # System-wide depth distribution (subsample for performance)
      if (nrow(values$raster_df) > 5000) {
        system_sample <- values$raster_df[sample(nrow(values$raster_df), 5000), ]
      } else {
        system_sample <- values$raster_df
      }

      # Sampled areas (areas with >5% DE)
      sampled_areas <- de_data[de_data$detection_efficiency > 0.05, ]

      # Create comparison dataframe for histogram
      system_depths <- data.frame(depth = system_sample$value, type = "System Total")
      sampled_depths <- data.frame(depth = sampled_areas$raster_value, type = "Effectively Sampled")
      comparison_data <- rbind(system_depths, sampled_depths)

      # Main histogram plot
      p1 <- ggplot(comparison_data, aes(x = depth, fill = type)) +
        geom_histogram(alpha = 0.7, position = "identity", bins = 30) +
        scale_fill_manual(values = c("System Total" = "lightblue", "Effectively Sampled" = "darkblue")) +
        labs(
          title = "Depth Distribution: System vs. Sampled Areas (>5% DE)",
          x = "Depth (m)",
          y = "Frequency",
          fill = "Coverage Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 9),
          legend.position = "bottom"
        )

      # Calculate differences for each depth bin
      depth_range <- range(c(system_depths$depth, sampled_depths$depth), na.rm = TRUE)
      breaks <- seq(depth_range[1], depth_range[2], length.out = 31)

      system_hist <- hist(system_depths$depth, breaks = breaks, plot = FALSE)
      sampled_hist <- hist(sampled_depths$depth, breaks = breaks, plot = FALSE)

      # Calculate proportional differences (sampled % - system %)
      system_props <- system_hist$counts / sum(system_hist$counts) * 100
      sampled_props <- sampled_hist$counts / sum(sampled_hist$counts) * 100
      diff_props <- sampled_props - system_props

      # Create difference dataframe
      bin_centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
      diff_data <- data.frame(
        depth = bin_centers,
        difference = diff_props,
        bias_type = ifelse(diff_props > 0, "Over-sampled", "Under-sampled")
      )

      # Difference plot
      p2 <- ggplot(diff_data, aes(x = depth, y = difference, fill = bias_type)) +
        geom_col(width = diff(breaks)[1] * 0.8) +
        geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
        scale_fill_manual(values = c("Over-sampled" = "blue", "Under-sampled" = "red")) +
        labs(
          title = "Sampling Bias by Depth",
          x = "Depth (m)",
          y = "Difference (%)",
          fill = "Bias Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 9),
          legend.position = "bottom"
        )

      # Return just the main histogram plot
      p1

    } else {
      # Placeholder plot
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Calculate distances and fit model\nto view depth coverage analysis",
                 size = 4, color = "gray50") +
        theme_void()
    }
  })

  # Sampling bias by depth plot
  output$sampling_bias_plot <- renderPlot({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data) && !is.null(values$raster_df)) {
      library(ggplot2)

      # System-wide depth distribution (subsample for performance)
      if (nrow(values$raster_df) > 5000) {
        system_sample <- values$raster_df[sample(nrow(values$raster_df), 5000), ]
      } else {
        system_sample <- values$raster_df
      }

      # Sampled areas (areas with >5% DE)
      sampled_areas <- de_data[de_data$detection_efficiency > 0.05, ]

      # Create comparison dataframe for histogram
      system_depths <- data.frame(depth = system_sample$value, type = "System Total")
      sampled_depths <- data.frame(depth = sampled_areas$raster_value, type = "Effectively Sampled")

      # Calculate differences for each depth bin
      depth_range <- range(c(system_depths$depth, sampled_depths$depth), na.rm = TRUE)
      breaks <- seq(depth_range[1], depth_range[2], length.out = 31)

      system_hist <- hist(system_depths$depth, breaks = breaks, plot = FALSE)
      sampled_hist <- hist(sampled_depths$depth, breaks = breaks, plot = FALSE)

      # Calculate proportional differences (sampled % - system %)
      system_props <- system_hist$counts / sum(system_hist$counts) * 100
      sampled_props <- sampled_hist$counts / sum(sampled_hist$counts) * 100
      diff_props <- sampled_props - system_props

      # Create difference dataframe
      bin_centers <- (breaks[-1] + breaks[-length(breaks)]) / 2
      diff_data <- data.frame(
        depth = bin_centers,
        difference = diff_props,
        bias_type = ifelse(diff_props > 0, "Over-sampled", "Under-sampled")
      )

      # Difference plot
      ggplot(diff_data, aes(x = depth, y = difference, fill = bias_type)) +
        geom_col(width = diff(breaks)[1] * 0.8) +
        geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
        scale_fill_manual(values = c("Over-sampled" = "blue", "Under-sampled" = "red")) +
        labs(
          title = "Sampling Bias by Depth",
          x = "Depth (m)",
          y = "Difference (%)",
          fill = "Bias Type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 10),
          legend.position = "bottom"
        )

    } else {
      # Placeholder plot
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Calculate distances and fit model\nto view sampling bias analysis",
                 size = 4, color = "gray50") +
        theme_void()
    }
  })

  # Performance metrics table
  output$performance_metrics <- renderTable({
    de_data <- get_detection_efficiency_data()
    if (!is.null(de_data) && !is.null(values$station_distances)) {
      # Calculate various performance metrics
      de_stats <- summary(de_data$detection_efficiency)

      # Distance to nearest receiver statistics
      min_distances <- values$station_distances %>%
        group_by(x, y) %>%
        summarise(min_dist = min(cost_distance, na.rm = TRUE), .groups = 'drop') %>%
        pull(min_dist)

      dist_stats <- summary(min_distances[!is.infinite(min_distances)])

      metrics <- data.frame(
        Metric = c(
          "DE Mean", "DE Median", "DE Std Dev",
          "Min Distance to Receiver (m)", "Mean Distance to Receiver (m)",
          "Number of Receivers", "System Coverage (%)"
        ),
        Value = c(
          round(de_stats["Mean"], 3),
          round(de_stats["Median"], 3),
          round(sd(de_data$detection_efficiency, na.rm = TRUE), 3),
          round(dist_stats["Min."], 0),
          round(dist_stats["Mean"], 0),
          nrow(values$receivers),
          round((sum(de_data$detection_efficiency > 0.05) / nrow(de_data)) * 100, 1)
        ),
        stringsAsFactors = FALSE
      )

      metrics
    } else {
      data.frame(
        Metric = c("DE Mean", "DE Median", "System Coverage (%)"),
        Value = c("--", "--", "--"),
        stringsAsFactors = FALSE
      )
    }
  }, striped = TRUE, hover = TRUE)
}
