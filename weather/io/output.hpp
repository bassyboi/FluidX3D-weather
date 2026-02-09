#pragma once
#include "../grid.hpp"
#include <string>

namespace io {
// Write VTK structured grid file for visualization in ParaView
void write_vtk(const Grid& grid, const std::string& dir, int step);

// Write CSV summary of current state (surface and column data)
void write_csv_summary(const Grid& grid, const std::string& dir, int step);

// Ensure output directory exists
void ensure_directory(const std::string& dir);
}
