# Load required packages
library(ismev)     # Introduction to Statistical Modeling of Extreme Values
library(car)       # For deltaMethod
library(dplyr)
library(readr)
library(tidyr)

# Load data
data <- read_csv("MAM_data.csv", show_col_types = FALSE)

# Convert to long format - handle missing values properly
data_long <- data %>%
  pivot_longer(-station, names_to = "year", values_to = "value") %>%
  mutate(year = as.numeric(year)) %>%
  filter(!is.na(value))  # Remove rows with missing values

# Function to process one station using ismev
process_station <- function(station_id, df) {
  cat("Processing station:", station_id, "\n")
  
  values <- df$value[df$station == station_id]
  values <- values[!is.na(values) & is.finite(values)]  # Remove NA and infinite values
  
  if (length(values) < 10) {
    warning("Insufficient data for station ", station_id, " (only ", length(values), " valid observations)")
    return(data.frame(
      station = station_id,
      location_est = NA, location_lower = NA, location_upper = NA, location_se = NA,
      scale_est = NA, scale_lower = NA, scale_upper = NA, scale_se = NA,
      shape_est = NA, shape_lower = NA, shape_upper = NA, shape_se = NA,
      nllh = NA,
      method = "insufficient_data"
    ))
  }
  
  tryCatch({
    # Fit GEV using ismev
    fit <- ismev::gev.fit(values, show = FALSE)
    
    # Extract parameters
    params <- fit$mle
    names(params) <- c("location", "scale", "shape")
    
    # Check for valid parameters
    if (any(is.na(params)) || any(!is.finite(params)) || params[2] <= 0) {
      warning("Invalid parameters for station ", station_id)
      return(data.frame(
        station = station_id,
        location_est = params[1], location_lower = NA, location_upper = NA, location_se = NA,
        scale_est = params[2], scale_lower = NA, scale_upper = NA, scale_se = NA,
        shape_est = params[3], shape_lower = NA, shape_upper = NA, shape_se = NA,
        nllh = fit$nllh,
        method = "invalid_parameters"
      ))
    }
    
    # Extract variance-covariance matrix
    vcov_matrix <- fit$cov
    rownames(vcov_matrix) <- colnames(vcov_matrix) <- names(params)
    
    # Check if variance-covariance matrix is valid
    if (is.null(vcov_matrix) || any(is.na(vcov_matrix)) || any(!is.finite(vcov_matrix))) {
      warning("Invalid variance-covariance matrix for station ", station_id)
      return(data.frame(
        station = station_id,
        location_est = params[1], location_lower = NA, location_upper = NA, location_se = NA,
        scale_est = params[2], scale_lower = NA, scale_upper = NA, scale_se = NA,
        shape_est = params[3], shape_lower = NA, shape_upper = NA, shape_se = NA,
        nllh = fit$nllh,
        method = "ismev_no_ci"
      ))
    }
    
    # Check if matrix is positive definite
    eigenvals <- eigen(vcov_matrix, only.values = TRUE)$values
    if (any(eigenvals <= 1e-12)) {
      warning("Variance-covariance matrix is not positive definite for station ", station_id)
      return(data.frame(
        station = station_id,
        location_est = params[1], location_lower = NA, location_upper = NA, location_se = NA,
        scale_est = params[2], scale_lower = NA, scale_upper = NA, scale_se = NA,
        shape_est = params[3], shape_lower = NA, shape_upper = NA, shape_se = NA,
        nllh = fit$nllh,
        method = "ismev_singular_matrix"
      ))
    }
    
    # Compute confidence intervals using car::deltaMethod
    ci_location <- car::deltaMethod(params, "location", vcov_matrix, level = 0.95)
    ci_scale <- car::deltaMethod(params, "scale", vcov_matrix, level = 0.95)
    ci_shape <- car::deltaMethod(params, "shape", vcov_matrix, level = 0.95)
    print(ci_shape)
    cat("  Success with ismev + car::deltaMethod\n")
    
    # Return results
    data.frame(
      station = station_id,
      location_est = params[1],
      location_lower = ci_location[1, 3],
      location_upper = ci_location[1, 4],
      location_se = ci_location[1, 2],
      scale_est = params[2],
      scale_lower = ci_scale[1, 3],
      scale_upper = ci_scale[1, 4],
      scale_se = ci_scale[1, 2],
      shape_est = params[3],
      shape_lower = ci_shape[1, 3],
      shape_upper = ci_shape[1, 4],
      shape_se = ci_shape[1, 2],
      nllh = fit$nllh,
      method = "ismev_car_delta"
    )
    
  }, error = function(e) {
    warning("Failed to process station ", station_id, ": ", e$message)
    return(data.frame(
      station = station_id,
      location_est = NA, location_lower = NA, location_upper = NA, location_se = NA,
      scale_est = NA, scale_lower = NA, scale_upper = NA, scale_se = NA,
      shape_est = NA, shape_lower = NA, shape_upper = NA, shape_se = NA,
      nllh = NA,
      method = "failed"
    ))
  })
}

# Process all stations
stations <- unique(data_long$station)
cat("Processing", length(stations), "stations using ismev package...\n\n")

# Use base R instead of purrr::map_dfr
result_list <- lapply(stations, function(x) process_station(x, data_long))
results <- bind_rows(result_list[!sapply(result_list, is.null)])

# Display summary
cat("\nSummary of results:\n")
if(nrow(results) > 0) {
  print(table(results$method))
}

# Display results
cat("\nResults:\n")
print(results)

# Save results
write_csv(results, "MAM_Parm_ML.csv")
cat("\nResults saved to 'gev_parameters_ismev.csv'\n")

# Validation if true parameters exist
if(file.exists("true_gev_parameters.csv")) {
  true_params <- read_csv("true_gev_parameters.csv", show_col_types = FALSE)
  
  validation <- results %>%
    left_join(true_params, by = "station") %>%
    filter(!is.na(location_est)) %>%  # Only validate successful fits
    mutate(
      loc_bias = location_est - location,
      scale_bias = scale_est - scale,
      shape_bias = shape_est - shape,
      loc_coverage = ifelse(!is.na(location_lower) & !is.na(location_upper),
                            location >= location_lower & location <= location_upper, NA),
      scale_coverage = ifelse(!is.na(scale_lower) & !is.na(scale_upper),
                              scale >= scale_lower & scale <= scale_upper, NA),
      shape_coverage = ifelse(!is.na(shape_lower) & !is.na(shape_upper),
                              shape >= shape_lower & shape <= shape_upper, NA)
    )
  
  if(nrow(validation) > 0) {
    cat("\nValidation Results:\n")
    cat("Location - Mean Bias:", round(mean(validation$loc_bias, na.rm = TRUE), 4),
        "| RMSE:", round(sqrt(mean(validation$loc_bias^2, na.rm = TRUE)), 4),
        "| Coverage:", round(mean(validation$loc_coverage, na.rm = TRUE), 3), "\n")
    cat("Scale    - Mean Bias:", round(mean(validation$scale_bias, na.rm = TRUE), 4),
        "| RMSE:", round(sqrt(mean(validation$scale_bias^2, na.rm = TRUE)), 4),
        "| Coverage:", round(mean(validation$scale_coverage, na.rm = TRUE), 3), "\n")
    cat("Shape    - Mean Bias:", round(mean(validation$shape_bias, na.rm = TRUE), 4),
        "| RMSE:", round(sqrt(mean(validation$shape_bias^2, na.rm = TRUE)), 4),
        "| Coverage:", round(mean(validation$shape_coverage, na.rm = TRUE), 3), "\n")
    
    write_csv(validation, "gev_validation_ismev.csv")
    cat("\nValidation results saved to 'gev_validation_ismev.csv'\n")
    
    # Display parameter comparison
    cat("\nParameter Comparison (Estimated vs True):\n")
    comparison_summary <- validation %>%
      select(station, location_est, location, scale_est, scale, shape_est, shape) %>%
      rename(loc_est = location_est, loc_true = location,
             scale_est = scale_est, scale_true = scale,
             shape_est = shape_est, shape_true = shape)
    print(comparison_summary)
  } else {
    cat("\nNo successful fits to validate.\n")
  }
}