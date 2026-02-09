#pragma once
#include "../grid.hpp"
#include <string>

namespace io {
// Load atmospheric initial conditions from a simplified GFS CSV format
// Format: level_hPa, temperature_K, u_wind_ms, v_wind_ms, humidity_kgkg
void load_gfs(Grid& grid, const std::string& path);
}
