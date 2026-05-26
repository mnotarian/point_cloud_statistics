# =============================================================================
# Script accompanying Motz and Notarian "3D modeling of archaeological excavations 
# using mobile devices as a complement to high-quality photogrammetry." 
# Submitted to Advances in Archaeological Practice in May 2026
# -----------------------------------------------------------------------------
# Purpose : Computes summary statistics and generates histograms for point-cloud
#           comparison data exported from CloudCompare. Uses two types of measurements:
#             - C2C (Cloud-to-Cloud) absolute distances
#             - C2M (Cloud-to-Mesh) signed distances
#           Histograms are colored to match CloudCompare's scalar-field
#           rendering by averaging the per-point RGB values within each bin.
#
# Input   : Two CSV files exported from CloudCompare, one per metric. They
#           are loaded separately because we used different scalar-field
#           color ramps for each (C2C: blue-green-yellow-red-purple; C2M:
#           diverging blue-white-red scheme). In CloudCompare, the scalar
#           field colors must be exported to RGB before each point cloud is
#           saved as ASCII CSV.
#
# Script Author  : Matthew Notarian
# First version, used for analysis as presented in the article, created by author
# Prepared for publication (cleaned, annotated) in Claude Opus 4.7
# Date    : 5-21-25
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Dependencies
# -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)


# -----------------------------------------------------------------------------
# 2. User-configurable parameters
# -----------------------------------------------------------------------------
# Stratigraphic Unit Label / label used in plot titles and output filenames.
su_label    <- "SU 13010"

# Number of histogram bins. A high bin count is used so that the per-bin
# averaged colors reproduce the CloudCompare gradient smoothly.
num_bins    <- 1000

# Output directory and dimensions for the saved figure.
output_dir  <- "C:/path/to/your/folder"
plot_width  <- 15   # cm
plot_height <- 15   # cm

# Input file paths. Set to NULL to be prompted interactively, or assign a
# path string for reproducible runs (recommended for final figure generation
# accompanying publication).
c2c_file <- NULL   # e.g. "C:/your/path/SU_13010_C2C.csv"
c2m_file <- NULL   # e.g. "C:/your/path/SU_13010_C2M.csv"


# -----------------------------------------------------------------------------
# 3. Data loading
# -----------------------------------------------------------------------------

#' Load a CloudCompare CSV export, prompting interactively if no path given.
#'
#' @param path   File path, or NULL to prompt with rstudioapi::selectFile().
#' @param label  Short description shown in the file-picker caption.
#' @return       Data frame, or NULL if the user cancels the prompt.
load_cloudcompare_csv <- function(path = NULL, label = "file") {
  if (is.null(path) || !nzchar(path)) {
    path <- rstudioapi::selectFile(paste("Select", label))
  }
  if (is.null(path) || !nzchar(path)) return(NULL)
  read.csv(path)
}

# Each metric has its own CSV because each carries the RGB columns matched
# to its own CloudCompare color ramp.
c2c_data <- load_cloudcompare_csv(c2c_file, label = "C2C CSV")
c2m_data <- load_cloudcompare_csv(c2m_file, label = "C2M CSV")


# -----------------------------------------------------------------------------
# 4. Summary statistics
# -----------------------------------------------------------------------------

#' Compute summary statistics for C2C absolute distances.
#'
#' Reports central tendency (mean, SD, median) and the cumulative percentage
#' of points falling within thresholds
#' (1 cm, 2 cm, 5 cm, 10 cm), along with the percentage of points in each
#' inter-threshold band and beyond 10 cm.
#'
#' @param x Numeric vector of C2C absolute distances (in meters).
#' @return  Named list of summary statistics, or NULL if x has no non-NA values.
summarize_c2c <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    warning("No non-missing C2C values; skipping C2C summary.")
    return(NULL)
  }
  cdf <- ecdf(x)
  pct_within_1cm  <- 100 * cdf(0.01)
  pct_within_2cm  <- 100 * cdf(0.02)
  pct_within_5cm  <- 100 * cdf(0.05)
  pct_within_10cm <- 100 * cdf(0.10)

  list(
    mean             = mean(x),
    sd               = sd(x),
    median           = median(x),
    pct_within_1cm   = pct_within_1cm,
    pct_within_2cm   = pct_within_2cm,
    pct_1cm_to_2cm   = pct_within_2cm  - pct_within_1cm,
    pct_2cm_to_5cm   = pct_within_5cm  - pct_within_2cm,
    pct_5cm_to_10cm  = pct_within_10cm - pct_within_5cm,
    pct_beyond_10cm  = 100 - pct_within_10cm
  )
}

#' Compute summary statistics for C2M signed distances.
#'
#' Reports mean and SD along with the percentage of positive vs. negative
#' deviations (i.e., the proportion of points above and below the reference
#' mesh surface).
#'
#' @param x Numeric vector of C2M signed distances (in meters).
#' @return  Named list of summary statistics, or NULL if x has no non-NA values.
summarize_c2m <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    warning("No non-missing C2M values; skipping C2M summary.")
    return(NULL)
  }
  pct_negative <- 100 * ecdf(x)(0)
  list(
    mean         = mean(x),
    sd           = sd(x),
    pct_positive = 100 - pct_negative,
    pct_negative = pct_negative
  )
}


# -----------------------------------------------------------------------------
# 5. Histogram construction
# -----------------------------------------------------------------------------
# Histograms are colored by averaging the CloudCompare scalar-field RGB values
# within each bin, producing a fill gradient that matches the source rendering.

#' Build the plotting data frame and per-bin averaged hex colors.
#'
#' @param distances Numeric vector of distance values for the x-axis.
#' @param r,g,b     Numeric vectors of per-point RGB values (0-255 scale).
#' @param n_bins    Number of histogram bins.
#' @return          List with $plot_data (binned points joined to colors) and
#'                  $bin_colors (per-bin averaged hex colors, in bin order).
prepare_histogram_data <- function(distances, r, g, b, n_bins) {
  data <- data.frame(value = distances, R = r, G = g, B = b)

  bin_breaks <- seq(min(data$value, na.rm = TRUE),
                    max(data$value, na.rm = TRUE),
                    length.out = n_bins + 1)
  data <- data %>%
    mutate(bin = cut(value, breaks = bin_breaks, include.lowest = TRUE))

  bin_colors <- data %>%
    group_by(bin) %>%
    summarise(
      R_avg = mean(R) / 255,
      G_avg = mean(G) / 255,
      B_avg = mean(B) / 255,
      .groups = "drop"
    ) %>%
    mutate(hex_color = rgb(R_avg, G_avg, B_avg))

  list(
    plot_data  = data %>% left_join(bin_colors, by = "bin"),
    bin_colors = bin_colors
  )
}

#' Render a color-matched histogram.
#'
#' Produces an image histogram with major/minor tick marks at
#' 5 cm / 1 cm intervals, axis lines drawn through zero and the data minimum,
#' and bin fill colors averaged from the CloudCompare scalar field.
#'
#' @param hist_data List returned by prepare_histogram_data().
#' @param title     Plot title.
#' @param x_label   x-axis label.
#' @param signed    Logical. TRUE for C2M (signed) data, which extends below
#'                  zero and uses the data minimum as the left axis line.
#'                  FALSE for C2C (absolute) data, which is bounded at zero.
#' @return          A ggplot object.
build_histogram <- function(hist_data, title, x_label = "Distance (m)",
                            signed = FALSE) {
  plot_data  <- hist_data$plot_data
  bin_colors <- hist_data$bin_colors

  data_min <- min(plot_data$value, na.rm = TRUE)
  data_max <- max(plot_data$value, na.rm = TRUE)

  # Anchor major/minor tick positions to clean 0.05 m / 0.01 m intervals.
  major_breaks <- seq(floor(data_min   / 0.05) * 0.05,
                      ceiling(data_max / 0.05) * 0.05, by = 0.05)
  minor_breaks <- seq(floor(data_min   / 0.01) * 0.01,
                      ceiling(data_max / 0.01) * 0.01, by = 0.01)

  # For signed data, draw the left vertical axis at the data minimum;
  # for absolute (non-negative) data, draw it at zero.
  left_axis_x <- if (signed) data_min else 0

  ggplot(plot_data, aes(x = value)) +
    geom_histogram(aes(fill = bin), bins = length(bin_colors$hex_color)) +
    scale_fill_manual(values = bin_colors$hex_color, guide = "none") +
    geom_vline(xintercept = left_axis_x) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(
      breaks       = major_breaks,
      minor_breaks = minor_breaks,
      guide        = guide_axis(minor.ticks = TRUE),
      expand       = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title, x = x_label, y = "Count") +
    theme_minimal() +
    theme(
      plot.title         = element_text(hjust = 0.5),
      panel.grid.major   = element_line(linewidth = 0.8),
      panel.grid.minor   = element_line(linewidth = 0.4),
      axis.ticks         = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length  = unit(0.2, "cm")
    )
}


# -----------------------------------------------------------------------------
# 6. Generate and save outputs
# -----------------------------------------------------------------------------

# --- C2C absolute distances --------------------------------------------------
if (!is.null(c2c_data)) {
  c2c_stats <- summarize_c2c(c2c_data$C2C.absolute.distances.Quadric..k.6.)
  print(c2c_stats)

  c2c_hist_data <- prepare_histogram_data(
    distances = c2c_data$C2C.absolute.distances.Quadric..k.6.,
    r = c2c_data$R, g = c2c_data$G, b = c2c_data$B,
    n_bins = num_bins
  )
  c2c_plot <- build_histogram(
    c2c_hist_data,
    title  = paste(su_label, "- C2C Absolute Distances"),
    signed = FALSE
  )
  print(c2c_plot)

  ggsave(
    filename = paste0(su_label, "_C2C_histogram.png"),
    plot     = c2c_plot,
    path     = output_dir,
    width    = plot_width,
    height   = plot_height,
    units    = "cm",
    bg       = "white"
  )
}

# --- C2M signed distances ----------------------------------------------------
if (!is.null(c2m_data)) {
  c2m_stats <- summarize_c2m(c2m_data$C2M.signed.distances)
  print(c2m_stats)

  c2m_hist_data <- prepare_histogram_data(
    distances = c2m_data$C2M.signed.distances,
    r = c2m_data$R, g = c2m_data$G, b = c2m_data$B,
    n_bins = num_bins
  )
  c2m_plot <- build_histogram(
    c2m_hist_data,
    title  = paste(su_label, "- C2M Signed Distances"),
    signed = TRUE
  )
  print(c2m_plot)

  ggsave(
    filename = paste0(su_label, "_C2M_histogram.png"),
    plot     = c2m_plot,
    path     = output_dir,
    width    = plot_width,
    height   = plot_height,
    units    = "cm",
    bg       = "white"
  )
}

