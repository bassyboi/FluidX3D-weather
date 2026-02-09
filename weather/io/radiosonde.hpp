#pragma once
#include "../grid.hpp"
#include <string>

namespace io {
// Load radiosonde vertical profile data
// Format CSV: pressure_hPa, height_m, temperature_C, dewpoint_C, wind_dir_deg, wind_speed_knots
void ingest_radiosonde(Grid& grid, const std::string& path);
}
