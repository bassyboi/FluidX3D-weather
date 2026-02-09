#pragma once
#include "../grid.hpp"
#include <string>

namespace io {
// Load surface weather station data
// Format CSV: time_s, temperature_C, humidity_pct, pressure_hPa, wind_speed_ms, wind_dir_deg
void ingest_aws(Grid& grid, const std::string& path);
}
