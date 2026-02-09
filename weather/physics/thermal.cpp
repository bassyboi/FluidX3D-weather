#include "thermal.hpp"
#include "../constants.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

namespace physics {

// ---------------------------------------------------------------------------
// Buoyancy force using Boussinesq approximation
// dw += g * (theta'/theta_ref + 0.608*qv - qc - qr)
// theta_ref is the horizontally-averaged theta at each vertical level
// ---------------------------------------------------------------------------
void apply_thermal_buoyancy(Grid& grid) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;

    // Compute horizontally-averaged theta at each vertical level (theta_ref)
    std::vector<double> theta_ref(Nz, 0.0);
    std::vector<int>    count(Nz, 0);

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (!grid.solid[id]) {
                    theta_ref[k] += grid.theta[id];
                    count[k]++;
                }
            }
        }
        if (count[k] > 0) {
            theta_ref[k] /= count[k];
        } else {
            // Fallback: use standard atmosphere value for this level
            theta_ref[k] = atm::T0 * std::pow(atm::p0 / grid.p[grid.idx(0, 0, k)], atm::kappa);
        }
    }

    // Apply buoyancy to vertical momentum tendency
    for (int k = 0; k < Nz; ++k) {
        double tref = theta_ref[k];
        if (tref <= 0.0) continue; // safety guard

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;

                double theta_prime = grid.theta[id] - tref;
                double qv_val = std::max(grid.qv[id], 0.0);
                double qc_val = std::max(grid.qc[id], 0.0);
                double qr_val = std::max(grid.qr[id], 0.0);

                double B = atm::g * (theta_prime / tref + 0.608 * qv_val - qc_val - qr_val);
                grid.dw[id] += B;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Newtonian cooling: relax theta toward an equilibrium profile
// Equilibrium uses a tropospheric lapse rate of 6.5 K/km, converted to
// potential temperature via the Exner function.
// Timescale tau_rad = 43200 s (12 hours).
// ---------------------------------------------------------------------------
void apply_radiation(Grid& grid, double /*dt*/) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;

    constexpr double lapse_rate = 6.5e-3;   // K/m (tropospheric)
    constexpr double tau_rad    = 43200.0;   // relaxation timescale [s] (12 hours)

    for (int k = 0; k < Nz; ++k) {
        double z_k = grid.z[k];

        // Equilibrium temperature: simple tropospheric lapse rate
        // Clamp at a minimum tropopause temperature (~200 K)
        double T_eq = std::max(atm::T0 - lapse_rate * z_k, 200.0);

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;

                // Convert T_eq to potential temperature at the local pressure
                double p_local = grid.p[id];
                double theta_eq = T_eq * std::pow(atm::p0 / p_local, atm::kappa);

                // Newtonian relaxation: dtheta/dt = -(theta - theta_eq) / tau_rad
                grid.dtheta[id] += -(grid.theta[id] - theta_eq) / tau_rad;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Surface heat flux: applies shortwave-like heating to the lowest level (k=0)
// Q = solar_constant * max(0, cos(solar_angle)) * surface_absorption_factor
// dtheta += Q / (rho * cp_d * dz[0])
// ---------------------------------------------------------------------------
void surface_heat_flux(Grid& grid, double solar_angle_rad, double /*dt*/) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;

    // Net shortwave absorbed at surface (accounts for albedo, atmospheric absorption, etc.)
    constexpr double surface_absorption = 0.3;
    double cos_angle = std::cos(solar_angle_rad);
    double Q = atm::solar_constant * std::max(0.0, cos_angle) * surface_absorption;

    if (Q <= 0.0) return; // nighttime, no heating

    double dz0 = grid.dz[0];

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            std::size_t id = grid.idx(i, j, 0);
            if (grid.solid[id]) continue;

            double rho_local = grid.rho[id];
            if (rho_local <= 0.0) continue; // safety guard

            // Heat flux warms the lowest grid cell
            grid.dtheta[id] += Q / (rho_local * atm::cp_d * dz0);
        }
    }
}

} // namespace physics
