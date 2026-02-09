#pragma once
#include "grid.hpp"
#include <string>

// Load terrain heightmap from a PGM (P2/P5) image or generate procedural terrain
void load_terrain(Grid& grid, const std::string& path);

// Generate simple procedural terrain (Gaussian hills)
void generate_terrain(Grid& grid, double max_height);

// Apply terrain to grid: mark cells below terrain as solid
void apply_terrain_to_grid(Grid& grid);
