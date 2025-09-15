#' Sample points from system detection efficiency surfaces
#'
#' Generates random points based on cumulative probability surfaces from system
#' detection efficiency calculations. Points are sampled with probability proportional
#' to the cumulative detection probability across the receiver array.
#'
#' @param system_de A list returned by detection efficiency system calculation
#'   containing a data tibble with cumulative_prob values.
#' @param n_points Integer. Number of points to sample per fish-time combination
#'   (if fish_id and time columns exist) or total. Default is 100.
#' @param prob_column Character. Name of the probability column to use for sampling.
#'   Default is "cumulative_prob".
#' @param fish_select Integer, character vector, or NULL. Fish ID(s) to sample from.
#'   Only used if fish_id column exists in the data. If NULL, samples from all fish.
#'   Default is NULL.
#' @param time_select Numeric vector, character vector, POSIXct vector, or NULL.
#'   Time period(s) to sample from. Only used if time-related columns exist.
#'   If NULL, samples from all time periods. Default is NULL.
#' @param position_points An sf object (e.g., from sample_points_from_probabilities) 
#'   containing fish_id and time_period information to use as a template. When provided,
#'   this overrides fish_select and time_select parameters, sampling the same fish-time
#'   combinations present in position_points. Default is NULL.
#' @param min_prob_threshold Numeric. Minimum probability threshold (0-1).
#'   Only cells with probability above this threshold are eligible for sampling.
#'   Default is 0.001 to exclude zero-probability cells.
#' @param max_prob_threshold Numeric. Maximum probability threshold (0-1).
#'   Only cells with probability below this threshold are eligible for sampling.
#'   Default is 1.0 (no upper limit).
#' @param seed Integer. Random seed for reproducible sampling. Default is NULL.
#' @param by_group Logical. If TRUE and fish/time columns exist, samples n_points
#'   for each fish-time combination. If FALSE or no grouping columns exist,
#'   samples n_points total. Default is TRUE.
#' @param uniform Logical. If TRUE, samples uniformly from all eligible cells.
#'   If FALSE, samples with probability proportional to prob_column values.
#'   Default is FALSE (probability-weighted sampling).
#' @param crs Coordinate reference system for the output sf object. Can be:
#'   \itemize{
#'     \item NULL (default) - attempts to use CRS from system_de$crs, falls back to WGS84
#'     \item Numeric EPSG code (e.g., 4326 for WGS84, 32618 for UTM Zone 18N)
#'     \item Character proj4 string
#'     \item An sf/sfc object from which to extract CRS
#'   }
#'
#' @return An sf object containing the sampled points with columns:
#'   \item{x}{X coordinates}
#'   \item{y}{Y coordinates}
#'   \item{probability}{The cumulative probability value used for sampling}
#'   \item{sample_id}{Sequential sample identifier}
#'   \item{fish_id}{Fish identifier (if present in input data)}
#'   \item{time_period}{Time period identifier (if present in input data)}
#'   \item{time_period_posix}{POSIXct datetime (if present in input data)}
#'   \item{time_period_label}{Human-readable time label (if present in input data)}
#'   \item{group_id}{Unique identifier for each fish-time combination (if applicable)}
#'   \item{geometry}{sf point geometry}
#'
#' @details
#' This function performs weighted random sampling where each spatial cell has
#' a probability of being selected proportional to its cumulative detection
#' probability value. This is useful for:
#' \itemize{
#'   \item Simulating animal release locations based on detection coverage
#'   \item Monte Carlo analysis of array performance
#'   \item Generating test positions weighted by detection probability
#'   \item Array design optimization studies
#' }
#'
#' The function automatically detects whether fish_id and time-related columns
#' exist in the input data and handles grouping accordingly. If these columns
#' are not present, it performs simple random sampling across all cells.
#'
#' @examples
#' \dontrun{
#' # Sample from system detection efficiency
#' sampled_points <- sample_points_from_system_de(
#'   system_DE,
#'   n_points = 500,
#'   seed = 123
#' )
#'
#' # Focus on moderate detection areas
#' moderate_samples <- sample_points_from_system_de(
#'   system_DE,
#'   n_points = 200,
#'   min_prob_threshold = 0.3,
#'   max_prob_threshold = 0.7
#' )
#'
#' # Sample with specific CRS (UTM Zone 18N)
#' utm_samples <- sample_points_from_system_de(
#'   system_DE,
#'   n_points = 300,
#'   crs = 32618  # EPSG code for UTM Zone 18N
#' )
#'
#' # Uniform sampling (ignore probability values)
#' uniform_samples <- sample_points_from_system_de(
#'   system_DE,
#'   n_points = 500,
#'   uniform = TRUE  # All eligible cells have equal probability
#' )
#'
#' # Plot sampled points
#' library(ggplot2)
#' ggplot() +
#'   geom_sf(data = sampled_points, aes(color = probability), alpha = 0.5) +
#'   scale_color_viridis_c(name = "Cumulative\nProbability") +
#'   theme_minimal()
#'
#' # If fish_id and time data exist in future versions:
#' # multi_samples <- sample_points_from_system_de(
#' #   system_DE,
#' #   n_points = 50,
#' #   fish_select = c(1, 2, 3),
#' #   time_select = c("2025-07-15", "2025-07-16"),
#' #   by_group = TRUE
#' # )
#' 
#' # Use position_points as a template for fish-time combinations
#' # First get position points from positioning results
#' position_points <- sample_points_from_probabilities(
#'   positioning_results,
#'   n_points = 100,
#'   fish_select = c(1, 2),
#'   time_select = c("2025-07-15", "2025-07-16")
#' )
#' 
#' # Then sample from system_de using same fish-time combinations
#' matched_de_samples <- sample_points_from_system_de(
#'   system_DE,
#'   n_points = 100,  # Same number per group as position_points
#'   position_points = position_points  # Use as template
#' )
#' }
#'
#' @seealso \code{\link{sample_points_from_probabilities}}
#'
#' @export
sample_points_from_system_de <- function(system_de,
                                        n_points = 100,
                                        prob_column = "cumulative_prob",
                                        fish_select = NULL,
                                        time_select = NULL,
                                        position_points = NULL,
                                        min_prob_threshold = 0.001,
                                        max_prob_threshold = 1.0,
                                        seed = NULL,
                                        by_group = TRUE,
                                        uniform = FALSE,
                                        crs = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Extract data from system_de object
  if ("data" %in% names(system_de)) {
    de_data <- as.data.frame(system_de$data)
  } else {
    stop("No 'data' element found in system_de object")
  }
  
  # Validate prob_column exists
  if (!prob_column %in% names(de_data)) {
    stop("Column '", prob_column, "' not found in system_de data. ",
         "Available columns: ", paste(names(de_data), collapse = ", "))
  }
  
  cat("Using probability column:", prob_column, "\n")
  cat("Original data:", nrow(de_data), "cells\n")
  
  # Handle position_points template if provided
  if (!is.null(position_points)) {
    if (!inherits(position_points, "sf")) {
      stop("position_points must be an sf object")
    }
    
    # Extract unique fish-time combinations from position_points
    template_data <- sf::st_drop_geometry(position_points)
    
    # Check what columns are available
    has_template_fish <- "fish_id" %in% names(template_data)
    has_template_time <- any(c("time_period", "time_period_posix", "time_period_label") %in% names(template_data))
    
    if (has_template_fish) {
      # Override fish_select with unique fish from template
      fish_select <- unique(template_data$fish_id)
      cat("Found", length(fish_select), "fish ID(s) in position_points template\n")
    }
    
    if (has_template_time) {
      # Determine which time column to use
      if ("time_period" %in% names(template_data)) {
        time_select <- unique(template_data$time_period)
        cat("Found", length(time_select), "time period(s) in position_points template\n")
      } else if ("time_period_label" %in% names(template_data)) {
        time_select <- unique(template_data$time_period_label)
        cat("Found", length(time_select), "time label(s) in position_points template\n")
      } else if ("time_period_posix" %in% names(template_data)) {
        time_select <- unique(template_data$time_period_posix)
        cat("Found", length(time_select), "time posix value(s) in position_points template\n")
      }
    }
    
    # Check if system_de data actually has these columns
    if ((has_template_fish && !("fish_id" %in% names(de_data))) ||
        (has_template_time && !any(c("time_period", "time_period_posix", "time_period_label") %in% names(de_data)))) {
      cat("Note: Template contains fish/time information, but system_de data lacks these columns.\n")
      if (by_group) {
        cat("      Will replicate sampling across template groups.\n")
      } else {
        cat("      Sampling will proceed without grouping.\n")
      }
    }
  }
  
  # Check for grouping columns
  has_fish_id <- "fish_id" %in% names(de_data)
  has_time <- any(c("time_period", "time_period_posix", "time_period_label") %in% names(de_data))
  
  # Filter by fish if specified and column exists
  if (!is.null(fish_select) && has_fish_id) {
    # Check if all fish_select values exist
    missing_fish <- setdiff(fish_select, unique(de_data$fish_id))
    if (length(missing_fish) > 0) {
      stop("Fish ID(s) not found: ", paste(missing_fish, collapse = ", "), 
           ". Available fish IDs: ", paste(unique(de_data$fish_id), collapse = ", "))
    }
    de_data <- de_data %>%
      dplyr::filter(fish_id %in% fish_select)
    cat("Filtered to fish:", paste(fish_select, collapse = ", "), 
        "-", nrow(de_data), "cells\n")
  } else if (!is.null(fish_select) && !has_fish_id) {
    if (!is.null(position_points)) {
      warning("position_points template contains fish_id, but system_de data has no fish_id column. Ignoring fish selection.")
    } else {
      warning("fish_select specified but no fish_id column found in data. Ignoring fish selection.")
    }
  }
  
  # Filter by time if specified and columns exist
  if (!is.null(time_select) && has_time) {
    # Handle different time selection formats
    if (is.numeric(time_select) && "time_period" %in% names(de_data)) {
      de_data <- de_data %>%
        dplyr::filter(time_period %in% time_select)
    } else if (is.character(time_select) && "time_period_label" %in% names(de_data)) {
      de_data <- de_data %>%
        dplyr::filter(time_period_label %in% time_select)
    } else if (inherits(time_select, "POSIXt") && "time_period_posix" %in% names(de_data)) {
      # Match by date
      time_chars <- as.character(as.Date(time_select))
      de_data <- de_data %>%
        dplyr::filter(as.character(as.Date(time_period_posix)) %in% time_chars)
    }
    cat("Filtered to time periods:", 
        ifelse("time_period" %in% names(de_data), 
               length(unique(de_data$time_period)), "unknown"), 
        "-", nrow(de_data), "cells\n")
  } else if (!is.null(time_select) && !has_time) {
    if (!is.null(position_points)) {
      warning("position_points template contains time information, but system_de data has no time columns. Ignoring time selection.")
    } else {
      warning("time_select specified but no time columns found in data. Ignoring time selection.")
    }
  }
  
  # Filter by probability thresholds
  de_data <- de_data %>%
    dplyr::filter(.data[[prob_column]] > min_prob_threshold, 
                  .data[[prob_column]] < max_prob_threshold,
                  !is.na(.data[[prob_column]]))
  
  cat("After probability thresholds (", min_prob_threshold, "< prob <", max_prob_threshold, "):", 
      nrow(de_data), "cells\n")
  
  if (nrow(de_data) == 0) {
    stop("No cells remain after filtering. Try adjusting min_prob_threshold (", min_prob_threshold, 
         ") or max_prob_threshold (", max_prob_threshold, ") or check your selection criteria.")
  }
  
  # Determine if we should process by groups
  # Process by groups if either:
  # 1. The data has fish/time columns, OR
  # 2. We have a position_points template with fish/time info (will replicate sampling)
  has_template_groups <- !is.null(position_points) && 
                        ((!is.null(fish_select) && length(fish_select) > 0) || 
                         (!is.null(time_select) && length(time_select) > 0))
  
  process_by_groups <- by_group && (has_fish_id || has_time || has_template_groups)
  
  if (process_by_groups && (has_fish_id || has_time)) {
    # Original logic for when data has grouping columns
    # Build grouping columns dynamically
    group_cols <- c()
    if (has_fish_id) group_cols <- c(group_cols, "fish_id")
    if ("time_period" %in% names(de_data)) group_cols <- c(group_cols, "time_period")
    if ("time_period_posix" %in% names(de_data)) group_cols <- c(group_cols, "time_period_posix")
    if ("time_period_label" %in% names(de_data)) group_cols <- c(group_cols, "time_period_label")
    
    # Get unique groups
    groups <- de_data %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))
    
    cat("\nProcessing", nrow(groups), "group combinations\n")
    
    # Initialize results list
    sampled_list <- list()
    
    # Process each group
    for (i in 1:nrow(groups)) {
      # Get current group values
      current_group <- groups[i, ]
      
      # Filter data for this group
      group_data <- de_data
      for (col in group_cols) {
        if (col %in% names(current_group)) {
          group_data <- group_data %>%
            dplyr::filter(.data[[col]] == current_group[[col]])
        }
      }
      
      if (nrow(group_data) == 0) next
      
      # Determine number of points to sample
      if (by_group) {
        group_n_points <- n_points
      } else {
        # Distribute points proportionally
        total_prob_mass <- sum(group_data[[prob_column]])
        overall_prob_mass <- sum(de_data[[prob_column]])
        group_n_points <- round(n_points * total_prob_mass / overall_prob_mass)
        group_n_points <- max(1, group_n_points)
      }
      
      # Check if we need replacement
      replacement <- group_n_points > nrow(group_data)
      
      # Prepare sampling weights
      if (uniform) {
        # Uniform sampling - all cells have equal probability
        sampling_weights <- NULL
      } else {
        # Probability-weighted sampling
        sampling_weights <- group_data[[prob_column]]
        sampling_weights <- sampling_weights / sum(sampling_weights)
      }
      
      # Sample indices
      sampled_indices <- sample(
        1:nrow(group_data), 
        size = group_n_points, 
        replace = replacement,
        prob = sampling_weights  # NULL for uniform sampling
      )
      
      # Extract sampled points
      sampled_data <- group_data[sampled_indices, ] %>%
        dplyr::select(x, y, probability = !!rlang::sym(prob_column),
                      dplyr::any_of(c("fish_id", "time_period", "time_period_posix", "time_period_label"))) %>%
        dplyr::mutate(
          sample_id = 1:group_n_points,
          group_id = paste(unlist(current_group), collapse = "_")
        )
      
      sampled_list[[i]] <- sampled_data
      
      # Progress message
      group_desc <- paste(names(current_group), "=", unlist(current_group), collapse = ", ")
      cat("  ", group_desc, ": sampled", group_n_points, "points from", nrow(group_data), "cells\n")
    }
    
    # Combine all sampled points
    all_sampled <- dplyr::bind_rows(sampled_list)
    
  } else if (process_by_groups && has_template_groups) {
    # Replicate sampling across template groups when data lacks grouping columns
    cat("\nReplicating sampling across template fish-time combinations\n")
    
    # Build template groups from position_points
    template_data <- sf::st_drop_geometry(position_points)
    
    # Get unique combinations based on what's available
    group_cols <- c()
    if (!is.null(fish_select)) group_cols <- c(group_cols, "fish_id")
    if (!is.null(time_select)) {
      if ("time_period" %in% names(template_data)) group_cols <- c(group_cols, "time_period")
      if ("time_period_posix" %in% names(template_data)) group_cols <- c(group_cols, "time_period_posix")
      if ("time_period_label" %in% names(template_data)) group_cols <- c(group_cols, "time_period_label")
    }
    
    # Get unique groups from template
    template_groups <- template_data %>%
      dplyr::select(dplyr::all_of(group_cols)) %>%
      dplyr::distinct() %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))
    
    cat("Found", nrow(template_groups), "unique fish-time combinations in template\n")
    
    # Initialize results list
    sampled_list <- list()
    
    # Sample once and replicate for each group
    # Determine sample size
    group_n_points <- n_points
    
    # Check if we need replacement
    replacement <- group_n_points > nrow(de_data)
    if (replacement) {
      warning("Requesting ", group_n_points, " points but only ", nrow(de_data), 
              " cells available. Sampling with replacement.")
    }
    
    # Prepare sampling weights
    if (uniform) {
      # Uniform sampling - all cells have equal probability
      sampling_weights <- NULL
    } else {
      # Probability-weighted sampling
      sampling_weights <- de_data[[prob_column]]
      sampling_weights <- sampling_weights / sum(sampling_weights)
    }
    
    # Sample indices (same for all groups)
    sampled_indices <- sample(
      1:nrow(de_data), 
      size = group_n_points, 
      replace = replacement,
      prob = sampling_weights  # NULL for uniform sampling
    )
    
    # Replicate for each template group
    for (i in 1:nrow(template_groups)) {
      # Get current group values
      current_group <- template_groups[i, ]
      
      # Extract sampled points
      sampled_data <- de_data[sampled_indices, ] %>%
        dplyr::select(x, y, probability = !!rlang::sym(prob_column))
      
      # Add group identifiers from template
      for (col in names(current_group)) {
        sampled_data[[col]] <- current_group[[col]]
      }
      
      # Add sample ID and group ID
      sampled_data <- sampled_data %>%
        dplyr::mutate(
          sample_id = 1:group_n_points,
          group_id = paste(unlist(current_group), collapse = "_")
        )
      
      sampled_list[[i]] <- sampled_data
      
      # Progress message
      group_desc <- paste(names(current_group), "=", unlist(current_group), collapse = ", ")
      cat("  ", group_desc, ": replicated", group_n_points, "points\n")
    }
    
    # Combine all sampled points
    all_sampled <- dplyr::bind_rows(sampled_list)
    
  } else {
    # Simple sampling without groups
    cat("\nProcessing all data as single group\n")
    
    # Check if we need replacement
    replacement <- n_points > nrow(de_data)
    
    if (replacement) {
      warning("Requesting ", n_points, " points but only ", nrow(de_data), 
              " cells available. Sampling with replacement.")
    }
    
    # Prepare sampling weights
    if (uniform) {
      # Uniform sampling - all cells have equal probability
      sampling_weights <- NULL
    } else {
      # Probability-weighted sampling
      sampling_weights <- de_data[[prob_column]]
      sampling_weights <- sampling_weights / sum(sampling_weights)
    }
    
    # Sample indices
    sampled_indices <- sample(
      1:nrow(de_data), 
      size = n_points, 
      replace = replacement,
      prob = sampling_weights  # NULL for uniform sampling
    )
    
    # Extract sampled points
    all_sampled <- de_data[sampled_indices, ] %>%
      dplyr::select(x, y, probability = !!rlang::sym(prob_column),
                    dplyr::any_of(c("fish_id", "time_period", "time_period_posix", "time_period_label"))) %>%
      dplyr::mutate(sample_id = 1:n_points)
    
    cat("Sampled", n_points, "points from", nrow(de_data), "cells\n")
  }
  
  # Re-number sample_id if needed
  if (!by_group && process_by_groups) {
    all_sampled$sample_id <- 1:nrow(all_sampled)
  }
  
  # Determine CRS to use
  if (!is.null(crs)) {
    # User specified CRS
    use_crs <- sf::st_crs(crs)
  } else if ("crs" %in% names(system_de)) {
    # Try to get from system_de
    use_crs <- sf::st_crs(system_de$crs)
  } else {
    # Default to WGS84
    use_crs <- sf::st_crs(4326)
    message("No CRS found in system_de object and none specified. Using WGS84 (EPSG:4326) as default.")
  }
  
  # Convert to sf object
  sampled_sf <- sf::st_as_sf(
    all_sampled, 
    coords = c("x", "y"), 
    crs = use_crs
  )
  
  # Add attributes for reference
  attr(sampled_sf, "prob_column") <- prob_column
  attr(sampled_sf, "n_points") <- n_points
  attr(sampled_sf, "min_prob_threshold") <- min_prob_threshold
  attr(sampled_sf, "max_prob_threshold") <- max_prob_threshold
  attr(sampled_sf, "by_group") <- by_group
  
  # Summary
  cat("\n=== SAMPLING SUMMARY ===\n")
  cat("Total points sampled:", nrow(sampled_sf), "\n")
  cat("Probability column:", prob_column, "\n")
  if (has_fish_id) {
    cat("Fish ID(s):", paste(unique(sampled_sf$fish_id), collapse = ", "), "\n")
  }
  if (has_time && "time_period_label" %in% names(sampled_sf)) {
    cat("Time periods:", length(unique(sampled_sf$time_period_label)), "\n")
  }
  if ("group_id" %in% names(sampled_sf)) {
    cat("Groups:", length(unique(sampled_sf$group_id)), "\n")
  }
  cat("Sampled probability range:", round(min(sampled_sf$probability), 4), 
      "to", round(max(sampled_sf$probability), 4), "\n")
  cat("Mean probability:", round(mean(sampled_sf$probability), 4), "\n")
  
  return(sampled_sf)
}