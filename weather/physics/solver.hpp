#pragma once
#include "../grid.hpp"

namespace physics {

// Advection using first-order upwind scheme for all prognostic variables
void advect(Grid& grid, double dt);

// Horizontal pressure gradient force: -1/rho * grad_h(p)
// Vertical PGF is hydrostatic; vertical motion comes from buoyancy + continuity
void pressure_gradient(Grid& grid);

// Coriolis force: f*(v, -u) where f = 2*Omega*sin(lat)
void coriolis(Grid& grid, double latitude_rad);

// Horizontal and vertical diffusion
void diffusion(Grid& grid, double Km, double Kh);

// Apply tendencies to update prognostic variables with forward Euler
void apply_tendencies(Grid& grid, double dt);

// Rayleigh damping sponge layer near domain top
void sponge_layer(Grid& grid, int sponge_depth);

// Enforce anelastic mass continuity: diagnose w from horizontal divergence
// w(k) = w(k-1) - (du/dx + dv/dy) * dz, with w(0) = 0 at surface
void enforce_continuity(Grid& grid);

// Clamp velocities to physically reasonable bounds for numerical safety
void clamp_velocities(Grid& grid, double max_horizontal, double max_vertical);

} // namespace physics
