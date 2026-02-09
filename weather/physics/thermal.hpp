#pragma once
#include "../grid.hpp"

namespace physics {

// Buoyancy force: vertical acceleration from temperature perturbation
// Uses Boussinesq approximation: dw += g * (theta'/theta_ref + 0.608*qv' - qc - qr)
void apply_thermal_buoyancy(Grid& grid);

// Simple radiation: Newtonian cooling toward equilibrium temperature profile
void apply_radiation(Grid& grid, double dt);

// Surface heat flux: warms lowest level based on solar angle
void surface_heat_flux(Grid& grid, double solar_angle_rad, double dt);

} // namespace physics
