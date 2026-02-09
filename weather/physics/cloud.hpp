#pragma once
#include "../grid.hpp"

namespace physics {

// Kessler warm-rain microphysics scheme
// Handles: condensation/evaporation, autoconversion, accretion, rain evaporation, sedimentation
void update_cloud_microphysics(Grid& grid, double dt);

} // namespace physics
