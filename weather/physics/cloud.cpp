#include "cloud.hpp"
#include "../constants.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

namespace physics {

// ---------------------------------------------------------------------------
// Kessler warm-rain microphysics scheme
//
// Steps applied per cell (non-solid only):
//   1. Saturation adjustment  (condensation / cloud evaporation)
//   2. Autoconversion         (cloud water -> rain water)
//   3. Accretion              (rain collecting cloud water)
//   4. Rain evaporation       (rain evaporating in sub-saturated air)
//   5. Rain sedimentation     (rain falling to lower cells)
//
// Saturation adjustment and latent heat effects modify grid fields directly
// (qv, qc, qr, theta) rather than through tendencies, because the
// condensation/evaporation step is a fast adjustment, not a slow tendency.
// ---------------------------------------------------------------------------

void update_cloud_microphysics(Grid& grid, double dt) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;

    // Temporary array to accumulate sedimentation flux into cells below.
    // sed_flux[idx] = rain mass flux arriving into cell idx from above [kg/kg].
    std::vector<double> sed_flux(grid.total(), 0.0);

    // Process from top to bottom so that sedimentation flux propagates downward
    for (int k = Nz - 1; k >= 0; --k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;

                double& qv_cell    = grid.qv[id];
                double& qc_cell    = grid.qc[id];
                double& qr_cell    = grid.qr[id];
                double& theta_cell = grid.theta[id];
                double  p_cell     = grid.p[id];
                double  rho_cell   = grid.rho[id];

                // Receive sedimentation flux from the cell above
                qr_cell += sed_flux[id];
                sed_flux[id] = 0.0; // consumed

                // Exner function: pi = (p/p0)^kappa, used for latent heat conversion
                double exner = std::pow(p_cell / atm::p0, atm::kappa);

                // Current temperature (needed for saturation calculation)
                double T_cell = theta_cell * exner;

                // -------------------------------------------------------
                // 1. Saturation adjustment (condensation / cloud evap)
                // -------------------------------------------------------
                double qvs = atm::saturation_mixing_ratio(p_cell, T_cell);

                if (qv_cell > qvs) {
                    // Supersaturated: condense excess vapor into cloud water
                    double delta_qc = qv_cell - qvs;
                    qv_cell -= delta_qc;
                    qc_cell += delta_qc;
                    // Latent heating: dtheta = Lv / (cp_d * exner) * delta_qc
                    theta_cell += (atm::Lv / (atm::cp_d * exner)) * delta_qc;
                } else if (qv_cell < qvs && qc_cell > 0.0) {
                    // Subsaturated with cloud present: evaporate cloud water
                    double deficit = qvs - qv_cell;
                    double delta_qc = std::min(qc_cell, deficit);
                    qv_cell += delta_qc;
                    qc_cell -= delta_qc;
                    // Latent cooling
                    theta_cell -= (atm::Lv / (atm::cp_d * exner)) * delta_qc;
                }

                // Update T after latent heat change for subsequent steps
                T_cell = theta_cell * exner;
                qvs = atm::saturation_mixing_ratio(p_cell, T_cell);

                // -------------------------------------------------------
                // 2. Autoconversion (cloud -> rain)
                //    Active when qc exceeds a threshold
                // -------------------------------------------------------
                if (qc_cell > atm::qc_autoconv_threshold) {
                    double dqr = atm::autoconv_rate * (qc_cell - atm::qc_autoconv_threshold) * dt;
                    dqr = std::min(dqr, qc_cell); // cannot convert more than available
                    qc_cell -= dqr;
                    qr_cell += dqr;
                }

                // -------------------------------------------------------
                // 3. Accretion (rain collecting cloud water)
                //    dqr = accretion_rate * qc * qr^0.875 * dt
                // -------------------------------------------------------
                if (qc_cell > 0.0 && qr_cell > 0.0) {
                    double dqr = atm::accretion_rate * qc_cell * std::pow(qr_cell, 0.875) * dt;
                    dqr = std::min(dqr, qc_cell);
                    qc_cell -= dqr;
                    qr_cell += dqr;
                }

                // -------------------------------------------------------
                // 4. Rain evaporation (in sub-saturated air below cloud)
                // -------------------------------------------------------
                if (qr_cell > 0.0 && qv_cell < qvs) {
                    double evap = atm::evaporation_rate * (qvs - qv_cell)
                                  * std::pow(qr_cell, 0.525) * dt;
                    evap = std::min(evap, qr_cell);
                    qr_cell -= evap;
                    qv_cell += evap;
                    // Latent cooling from rain evaporation
                    theta_cell -= (atm::Lv / (atm::cp_d * exner)) * evap;
                }

                // -------------------------------------------------------
                // 5. Sedimentation (rain falling)
                //    Vt = terminal_velocity_coeff * sqrt(rho0/rho) * qr^0.1375
                //    Fraction of rain that leaves this cell downward per dt
                // -------------------------------------------------------
                if (qr_cell > 0.0) {
                    double Vt = atm::terminal_velocity_coeff
                                * std::sqrt(atm::rho0 / std::max(rho_cell, 0.1))
                                * std::pow(qr_cell, 0.1375);

                    double dz_k = grid.dz[k];
                    // Fraction of rain leaving cell downward (CFL-limited to 1.0)
                    double frac = std::min(Vt * dt / dz_k, 1.0);
                    double qr_out = qr_cell * frac;

                    qr_cell -= qr_out;

                    if (k > 0) {
                        // Transfer to cell below via sedimentation flux array
                        std::size_t id_below = grid.idx(i, j, k - 1);
                        sed_flux[id_below] += qr_out;
                    }
                    // If k == 0, rain reaches the ground (accumulated as precipitation)
                    // Could be stored in a surface precipitation accumulator if desired
                }

                // -------------------------------------------------------
                // Clamp all moisture fields to non-negative values
                // -------------------------------------------------------
                qv_cell = std::max(qv_cell, 0.0);
                qc_cell = std::max(qc_cell, 0.0);
                qr_cell = std::max(qr_cell, 0.0);
            }
        }
    }
}

} // namespace physics
