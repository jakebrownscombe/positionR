library(shinydashboard)
library(shiny)
library(plotly)
library(DT)
library(ggplot2)
library(sf)

ui <- dashboardPage(
  dashboardHeader(title = "Array DesignR"),

  dashboardSidebar(
    width = 280,  # Smaller than 350px but bigger than default 230px
    sidebarMenu(
      menuItem("Array Design", tabName = "design", icon = icon("map"), selected = TRUE),
      menuItem("Range Modeling", tabName = "range", icon = icon("signal")),
      menuItem("Array Performance", tabName = "performance", icon = icon("chart-bar")),
      menuItem("Instructions", tabName = "instructions", icon = icon("info-circle"))
    ),

    hr(),

    # File input for raster
    fileInput("raster_file",
              "Load Raster File",
              accept = c(".tif", ".tiff", ".grd", ".asc")),

    hr(),

    # Regular point generation
    h4("Regular Array Generation"),

    numericInput("n_regular_points",
                 "Number of Points:",
                 value = 50,
                 min = 1,
                 max = 100,
                 step = 1),

    numericInput("array_seed",
                 "Random Seed:",
                 value = 123,
                 min = 1,
                 max = 999,
                 step = 1),

    actionButton("generate_regular",
                 "Generate",
                 icon = icon("th"),
                 class = "btn-info btn-sm",
                 style = "width: 90%; margin-bottom: 10px; font-size: 11px; padding: 4px 8px;"),

    actionButton("clear_receivers",
                 "Clear All",
                 icon = icon("trash"),
                 class = "btn-warning btn-sm",
                 style = "width: 90%; margin-bottom: 10px; font-size: 11px; padding: 4px 8px;"),

    actionButton("download_receivers",
                 "Download",
                 icon = icon("download"),
                 class = "btn-success btn-sm",
                 style = "width: 90%; margin-bottom: 10px; font-size: 11px; padding: 4px 8px;"),

    hr(),

    # Distance calculation
    h4("Distance Analysis"),

    numericInput("max_distance",
                 "Max Distance (m):",
                 value = 30000,
                 min = 1000,
                 max = 100000,
                 step = 1000),

    actionButton("calculate_distances",
                 "Calculate Distances",
                 icon = icon("calculator"),
                 class = "btn-info btn-sm",
                 style = "width: 90%; margin-bottom: 10px; font-size: 11px; padding: 4px 8px;"),

    checkboxInput("show_distances",
                  "Show Distance Fields",
                  value = FALSE),

    checkboxInput("show_detection_efficiency",
                  "Show Detection Efficiency",
                  value = FALSE),

    hr(),

    # Quick reference
    h4("Quick Reference"),
    tags$div(
      style = "font-size: 11px;",
      tags$ul(
        tags$li("Load raster → Place receivers → Calculate distances"),
        tags$li("Red dots: Manual receivers | Orange dots: Regular array"),
        tags$li("Toggle: Depth | Distance | Detection Efficiency views"),
        tags$li("See Instructions tab for detailed guidance")
      )
    )
  ),

  dashboardBody(
    tabItems(
      # Array Design tab
      tabItem(tabName = "design",
        fluidRow(
          box(
            title = "Interactive Array Design Map",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            height = "600px",

            plotlyOutput("raster_plot", height = "550px")
          )
        ),

        fluidRow(
          box(
            title = "Map Information",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            height = "250px",

            verbatimTextOutput("raster_info")
          ),

          box(
            title = "Current Receivers",
            status = "success",
            solidHeader = TRUE,
            width = 6,

            tableOutput("receiver_summary")
          )
        ),

        fluidRow(
          box(
            title = "Receiver Details",
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            DT::dataTableOutput("receiver_details")
          )
        )
      ),

      # Range Modeling tab
      tabItem(tabName = "range",
        fluidRow(
          box(
            title = "Detection Range Model Parameters",
            status = "primary",
            solidHeader = TRUE,
            width = 6,

            h4("Depth Range"),
            fluidRow(
              column(6, numericInput("min_depth", "Min Depth (m):", value = 1, min = 0.1, max = 50, step = 0.1)),
              column(6, numericInput("max_depth", "Max Depth (m):", value = 35, min = 1, max = 100, step = 1))
            ),

            h4("Detection Range at Shallow Depth"),
            fluidRow(
              column(6, numericInput("d50_min_depth", "50% Range (m):", value = 400, min = 50, max = 2000, step = 10)),
              column(6, numericInput("d95_min_depth", "5% Range (m):", value = 800, min = 100, max = 3000, step = 10))
            ),

            h4("Detection Range at Deep Depth"),
            fluidRow(
              column(6, numericInput("d50_max_depth", "50% Range (m):", value = 750, min = 50, max = 3000, step = 10)),
              column(6, numericInput("d95_max_depth", "5% Range (m):", value = 1500, min = 100, max = 5000, step = 10))
            ),

            br(),
            actionButton("fit_range_model", "Fit Detection Model",
                         icon = icon("chart-line"), class = "btn-primary btn-block")
          ),

          box(
            title = "Model Visualization",
            status = "info",
            solidHeader = TRUE,
            width = 6,

            plotOutput("range_model_plot", height = "400px")
          )
        ),

        fluidRow(
          box(
            title = "Detection Efficiency Predictions",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            height = "500px",

            plotlyOutput("detection_efficiency_plot", height = "450px")
          )
        )
      ),

      # Array Performance tab
      tabItem(tabName = "performance",
        # Row 1: Performance Metrics Summary
        fluidRow(
          valueBoxOutput("mean_de"),
          valueBoxOutput("coverage_area"),
          valueBoxOutput("high_de_area")
        ),

        # Row 2: Coverage Threshold Analysis
        fluidRow(
          box(
            title = "Detection Efficiency Distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            height = "400px",

            plotOutput("de_histogram", height = "350px")
          ),

          box(
            title = "Coverage Threshold & Performance Analysis",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            height = "400px",

            fluidRow(
              column(6,
                h5("Coverage Thresholds"),
                tableOutput("coverage_thresholds")
              ),
              column(6,
                h5("Array Performance Metrics"),
                tableOutput("performance_metrics")
              )
            )
          )
        ),

        # Row 3: Depth Coverage Analysis
        fluidRow(
          box(
            title = "Depth Coverage Comparison",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            height = "450px",

            plotOutput("depth_comparison", height = "400px")
          ),

          box(
            title = "Sampling Bias by Depth",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            height = "450px",

            plotOutput("sampling_bias_plot", height = "400px")
          )
        )
      ),

      # Instructions tab
      tabItem(tabName = "instructions",
        fluidRow(
          box(
            title = "Array DesignR - Scientific User Guide",
            status = "info",
            solidHeader = TRUE,
            width = 12,

            # Quick Start Section
            h3("🚀 Quick Start Guide"),
            tags$div(
              style = "background-color: #f8f9fa; padding: 15px; border-left: 4px solid #007bff; margin-bottom: 20px;",
              h4("Essential Workflow:"),
              tags$ol(
                tags$li(strong("Load raster data:"), " Upload .tif, .grd, or .asc bathymetry/habitat files"),
                tags$li(strong("Design array:"), " Place receivers manually (click map) or generate regular patterns"),
                tags$li(strong("Calculate distances:"), " Click 'Calculate Distances' to compute coverage"),
                tags$li(strong("Fit detection model:"), " Use Range Modeling tab (auto-fitted on startup)"),
                tags$li(strong("Assess performance:"), " Evaluate coverage and bias in Array Performance tab")
              )
            ),

            hr(),

            # Scientific Background Section
            h3("🔬 Scientific Background"),
            h4("Weighted Average Detection Efficiency (WADE) Methodology"),
            p("Array DesignR implements detection efficiency-based array design principles developed for acoustic telemetry systems.
              The approach integrates spatial detection probability modeling with cost-path distance calculations to optimize
              receiver placement for fish positioning studies."),

            h4("Key Concepts:"),
            tags$ul(
              tags$li(strong("Detection Efficiency (DE):"), " Probability of detecting a tagged animal at a given location,
                      incorporating distance-decay and environmental effects"),
              tags$li(strong("Cost-Path Distances:"), " Acoustic signal travel paths accounting for bathymetric barriers,
                      more realistic than Euclidean distances"),
              tags$li(strong("System-Wide Coverage:"), " Combined detection probability across all receivers,
                      calculated as 1 - ∏(1 - DE_i) for individual receiver efficiencies")
            ),

            hr(),

            # Module-Specific Instructions
            h3("📋 Module-Specific Guide"),

            # Array Design Module
            h4("🗺️ Array Design Module"),
            tags$div(
              style = "margin-left: 20px;",
              h5("Raster Data Requirements:"),
              tags$ul(
                tags$li("Formats: .tif (GeoTIFF), .grd (raster), .asc (ASCII grid)"),
                tags$li("Coordinate systems: UTM or geographic projections supported"),
                tags$li("Resolution: Optimize for study scale (typically 10-100m for lakes, 100-1000m for marine)")
              ),

              h5("Receiver Placement Strategies:"),
              tags$ul(
                tags$li(strong("Manual Placement:"), " Click map to add receivers at specific locations.
                        Optimal for irregular shorelines, known fish habitats, or targeted coverage"),
                tags$li(strong("Regular Arrays:"), " Generate systematic grids using positionR algorithms.
                        Optimal for broad coverage and unbiased sampling. Default: 50 receivers"),
                tags$li(strong("Hybrid Approach:"), " Combine regular arrays with manual adjustments for
                        optimal coverage while addressing specific study requirements")
              ),

              h5("Visualization Options:"),
              tags$ul(
                tags$li(strong("Depth View:"), " Default bathymetric visualization (blue gradient)"),
                tags$li(strong("Distance Fields:"), " Shows minimum distance to nearest receiver (red = close, yellow = far)"),
                tags$li(strong("Detection Efficiency:"), " System-wide detection probability surface using fitted model (magma colorscale)")
              )
            ),

            # Range Modeling Module
            h4("📊 Range Modeling Module"),
            tags$div(
              style = "margin-left: 20px;",
              h5("Detection Model Parameters:"),
              tags$ul(
                tags$li(strong("d50:"), " Distance at 50% detection probability - represents effective detection range"),
                tags$li(strong("d95:"), " Distance at 95% detection probability - represents minimum reliable detection range"),
                tags$li(strong("Depth Dependency:"), " Detection ranges typically increase with depth due to acoustic layering")
              ),

              h5("Parameter Considerations:"),
              tags$ul(
                tags$li("Base values on empirical range testing in similar environments"),
                tags$li("Shallow water (1-10m): Typically d50=200-400m, d95=400-800m"),
                tags$li("Deep water (>20m): Typically d50=500-1000m, d95=1000-2000m"),
                tags$li("Consider tag type, power settings, and environmental conditions")
              ),

              h5("Model Validation:"),
              p("Compare model predictions with actual range test data. Good models show realistic
                detection decay curves and align with known acoustic propagation physics.")
            ),

            # Performance Analysis Module
            h4("📈 Array Performance Module"),
            tags$div(
              style = "margin-left: 20px;",
              h5("Coverage Threshold Interpretation:"),
              tags$ul(
                tags$li(strong(">5% DE:"), " Minimum meaningful detection probability - defines 'monitored' area"),
                tags$li(strong(">25% DE:"), " Low but useable detection probability - suitable for presence/absence"),
                tags$li(strong(">50% DE:"), " Good detection probability - suitable for positioning and behavior"),
                tags$li(strong(">75% DE:"), " Excellent detection probability - high confidence in all metrics")
              ),

              h5("Sampling Bias Analysis:"),
              p("The depth distribution comparison identifies whether your array proportionally samples
                all available habitats. Red bars indicate under-sampled depth ranges that may bias
                ecological interpretations. Blue bars show over-sampled ranges."),

              h5("Performance Benchmarks:"),
              tags$ul(
                tags$li("Excellent arrays: >75% coverage at >5% DE, mean DE >0.4"),
                tags$li("Good arrays: >50% coverage at >5% DE, mean DE >0.25"),
                tags$li("Adequate arrays: >25% coverage at >5% DE, mean DE >0.15"),
                tags$li("Poor arrays: <25% coverage at >5% DE, significant depth bias")
              )
            ),

            hr(),

            # Best Practices Section
            h3("⚡ Best Practices & Considerations"),

            h4("Array Design Principles:"),
            tags$ul(
              tags$li("Receiver spacing should be 2-3× the effective detection range (d50)"),
              tags$li("Prioritize coverage of ecologically relevant habitats over uniform spacing"),
              tags$li("Consider acoustic shadowing from islands, shallow areas, or structure"),
              tags$li("Plan for equipment loss - include some redundancy in critical areas")
            ),

            h4("Study Integration:"),
            tags$ul(
              tags$li("Match array design to research questions (fine-scale behavior vs. broad movement)"),
              tags$li("Consider temporal aspects - seasonal habitat use, migration corridors"),
              tags$li("Integrate with other data collection (environmental sensors, fish sampling)"),
              tags$li("Plan for long-term maintenance and data retrieval")
            ),

            h4("Data Export & Analysis:"),
            tags$ul(
              tags$li("Export receiver coordinates as CSV for field deployment"),
              tags$li("Integrate with positionR package for WADE positioning analysis"),
              tags$li("Document array design decisions for publication methods"),
              tags$li("Validate array performance with field range testing")
            ),

            hr(),

            # Troubleshooting Section
            h3("🔧 Troubleshooting & Tips"),

            h4("Common Issues:"),
            tags$ul(
              tags$li(strong("Large raster files:"), " App automatically subsamples for performance.
                      Consider reducing resolution for initial design, then refine with higher resolution"),
              tags$li(strong("Model fitting errors:"), " Check parameter ranges are realistic.
                      d95 should be greater than d50, and both should reflect known acoustic ranges"),
              tags$li(strong("Slow distance calculations:"), " Large arrays (>50 receivers) or high-resolution
                      rasters may take several minutes. Progress is shown in notifications"),
              tags$li(strong("Memory issues:"), " Restart R session if calculations fail.
                      Consider reducing raster resolution or number of receivers")
            ),

            h4("Performance Optimization:"),
            tags$ul(
              tags$li("Start with moderate resolution rasters (1000×1000 pixels max)"),
              tags$li("Test array designs with ~20-30 receivers before scaling up"),
              tags$li("Use regular arrays as starting points, then manually adjust"),
              tags$li("Export and save intermediate results to avoid recalculation")
            ),

            hr(),

            # Integration Section
            h3("🔗 positionR Integration"),
            p("Array DesignR is designed to integrate seamlessly with the positionR package ecosystem:"),
            tags$ul(
              tags$li("Exported receiver coordinates can be imported directly into positionR workflows"),
              tags$li("Detection models fitted here are compatible with WADE positioning calculations"),
              tags$li("Performance metrics help validate positioning accuracy expectations"),
              tags$li("Array designs can be tested with simulated fish movements using positionR functions")
            ),

            tags$div(
              style = "background-color: #e8f5e8; padding: 15px; border-left: 4px solid #28a745; margin-top: 20px;",
              h4("📚 Further Reading:"),
              p("For detailed methodology and case studies, see the positionR package documentation and associated publications on WADE positioning theory.")
            )
          )
        )
      )
    )
  )
)