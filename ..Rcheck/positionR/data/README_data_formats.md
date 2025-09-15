# Data Format Guidelines

## Important: Excel Date Warning ⚠️

**CAUTION**: If you open the CSV files (`station_info.csv` or `temporal_info.csv`) in Excel, it will automatically convert the date formats and may corrupt the data when you save.

### The Problem:
- Original format: `2025-01-15` (YYYY-MM-DD)
- Excel converts to: `1/15/2025` or `15/1/2025` (depending on your region)
- This can cause confusion and data inconsistency

### Solutions:

#### Option 1: Use Text Import (Recommended)
1. Open Excel
2. Go to **Data > Get Data > From Text/CSV**
3. Select your CSV file
4. Choose **"Do not detect data types"** or set date columns to **Text**
5. This preserves the original YYYY-MM-DD format

#### Option 2: Function Handles Multiple Formats
The `simulate_fish_tracks()` function now automatically handles these common date formats:
- `2025-01-15` (ISO standard - preferred)
- `1/15/2025` (US Excel format)
- `15/1/2025` (European Excel format)
- `1-15-2025`, `15-1-2025` (alternative formats)
- `2025/01/15` (alternative ISO format)

#### Option 3: Edit Safely
If you must edit in Excel:
1. Make your changes
2. **Save As > CSV (Comma delimited)**
3. **DO NOT** use Excel's default save - it may change date formats

### Best Practice:
Keep the original YYYY-MM-DD format (`2025-01-15`) in your CSV files. The function will handle Excel-converted formats automatically, but the original format is clearest and most reliable.

## Data Structure

### species_movement_params.csv
Species-specific movement parameters for realistic fish movement simulation:
- `species`: Fish species name (Walleye, Smallmouth Bass, Muskellunge)
- `step_length_mean_base`: Base mean step length in meters
- `step_length_sd_base`: Base standard deviation of step length
- `turning_angle_mean`: Mean turning angle in degrees (typically 0)
- `turning_angle_sd`: Standard deviation of turning angle in degrees
- `size_scalar_slope`: Slope for size scaling (larger fish move more)
- `size_scalar_intercept`: Intercept for size scaling
- `min_size_cm`, `max_size_cm`, `typical_size_cm`: Size ranges in centimeters
- `movement_description`: Description of species movement behavior

**Size scaling formula:** `actual_step_length = base_step_length × (size_scalar_slope × fish_size_cm + size_scalar_intercept)`

### station_info.csv
- `station_id`: Unique identifier for each receiver station
- `x`, `y`: UTM coordinates of station location
- `depth_m`: Water depth at station location (meters)
- `start_date`: Date when station was deployed (YYYY-MM-DD)
- `end_date`: Date when station was retrieved (YYYY-MM-DD)

### temporal_info.csv
- `date`: Date for environmental conditions (YYYY-MM-DD)
- `water_temp_c`: Water temperature in Celsius
- `wind_speed_ms`: Wind speed in meters per second
- `turbidity_ntu`: Water turbidity in Nephelometric Turbidity Units
- `daylight_hours`: Hours of daylight
- `moon_phase`: Moon phase (0 = new moon, 1 = full moon)

## Usage Example

### Basic Species Simulation
```r
# Use species-specific movement parameters
fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances_df,
  species = "Walleye",           # Use Walleye movement parameters
  fish_size_cm = 50,             # 50cm fish (larger than typical = longer steps)
  n_paths = 10,
  n_steps = 100
)
```

### Manual Parameters (Traditional)
```r
# Specify individual movement parameters
fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances_df,
  step_length_mean = 60,         # Manual specification
  step_length_sd = 25,
  turning_angle_sd = 35
)
```

### Temporal DE with Species
```r
# Combine species parameters with temporal DE prediction
fish_simulation <- simulate_fish_tracks(
  raster = depth_raster,
  station_distances = station_distances_df,
  species = "Muskellunge",       # Large, directed movements
  fish_size_cm = 120,            # Large fish
  station_info = station_info,    # Excel dates OK
  temporal_info = temporal_info,  # Excel dates OK
  de_model = your_de_model
)
```

### Available Species
- **Walleye**: Moderate cruising with directed movements (45cm typical)
- **Smallmouth Bass**: Short bursts with frequent direction changes (35cm typical) 
- **Muskellunge**: Large directed movements with long glides (100cm typical)