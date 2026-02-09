#include "solver.hpp"
#include "../constants.hpp"
#include <cmath>
#include <algorithm>

namespace physics {

// ---------------------------------------------------------------------------
// Helper: first-order upwind advection for a single scalar field.
// Computes  -u * d(field)/dx  -  v * d(field)/dy  -  w * d(field)/dz
// and accumulates the result into the corresponding tendency array.
// ---------------------------------------------------------------------------
static void advect_field(Grid& grid,
                         const std::vector<double>& field,
                         std::vector<double>& tendency,
                         double /*dt*/) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;
    const double inv_dx = 1.0 / grid.dx;
    const double inv_dy = 1.0 / grid.dy;

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const std::size_t c = grid.idx(i, j, k);
                if (grid.solid[c]) continue;

                const double u_c = grid.u[c];
                const double v_c = grid.v[c];
                const double w_c = grid.w[c];
                const double f_c = field[c];

                // --- x-direction (periodic) ---
                double dfdx;
                if (u_c >= 0.0) {
                    // upwind from the left
                    const std::size_t im = grid.idx(grid.wrap_x(i - 1), j, k);
                    dfdx = (f_c - field[im]) * inv_dx;
                } else {
                    // upwind from the right
                    const std::size_t ip = grid.idx(grid.wrap_x(i + 1), j, k);
                    dfdx = (field[ip] - f_c) * inv_dx;
                }

                // --- y-direction (periodic) ---
                double dfdy;
                if (v_c >= 0.0) {
                    const std::size_t jm = grid.idx(i, grid.wrap_y(j - 1), k);
                    dfdy = (f_c - field[jm]) * inv_dy;
                } else {
                    const std::size_t jp = grid.idx(i, grid.wrap_y(j + 1), k);
                    dfdy = (field[jp] - f_c) * inv_dy;
                }

                // --- z-direction (clamped at boundaries) ---
                double dfdz;
                if (w_c >= 0.0) {
                    const int km = grid.clamp_z(k - 1);
                    if (km == k) {
                        dfdz = 0.0; // at bottom boundary
                    } else {
                        const double dz_below = 0.5 * (grid.dz[km] + grid.dz[k]);
                        dfdz = (f_c - field[grid.idx(i, j, km)]) / dz_below;
                    }
                } else {
                    const int kp = grid.clamp_z(k + 1);
                    if (kp == k) {
                        dfdz = 0.0; // at top boundary
                    } else {
                        const double dz_above = 0.5 * (grid.dz[k] + grid.dz[kp]);
                        dfdz = (field[grid.idx(i, j, kp)] - f_c) / dz_above;
                    }
                }

                tendency[c] += -(u_c * dfdx + v_c * dfdy + w_c * dfdz);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// advect: apply upwind advection to every prognostic variable
// ---------------------------------------------------------------------------
void advect(Grid& grid, double dt) {
    advect_field(grid, grid.u,     grid.du,     dt);
    advect_field(grid, grid.v,     grid.dv,     dt);
    advect_field(grid, grid.w,     grid.dw,     dt);
    advect_field(grid, grid.theta, grid.dtheta, dt);
    advect_field(grid, grid.qv,   grid.dqv,    dt);
    advect_field(grid, grid.qc,   grid.dqc,    dt);
    advect_field(grid, grid.qr,   grid.dqr,    dt);
}

// ---------------------------------------------------------------------------
// pressure_gradient: horizontal only (hydrostatic approximation)
//   du += -1/rho * dp/dx,  dv += -1/rho * dp/dy
// The vertical pressure gradient is in hydrostatic balance with gravity,
// so vertical accelerations come solely from the buoyancy perturbation
// (handled by apply_thermal_buoyancy). This is the standard approach
// for mesoscale NWP models at O(1km) resolution.
// ---------------------------------------------------------------------------
void pressure_gradient(Grid& grid) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;
    const double inv_2dx = 1.0 / (2.0 * grid.dx);
    const double inv_2dy = 1.0 / (2.0 * grid.dy);

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const std::size_t c = grid.idx(i, j, k);
                if (grid.solid[c]) continue;

                const double inv_rho = 1.0 / grid.rho[c];

                // x-gradient (periodic)
                const std::size_t ip = grid.idx(grid.wrap_x(i + 1), j, k);
                const std::size_t im = grid.idx(grid.wrap_x(i - 1), j, k);
                const double dpdx = (grid.p[ip] - grid.p[im]) * inv_2dx;

                // y-gradient (periodic)
                const std::size_t jp = grid.idx(i, grid.wrap_y(j + 1), k);
                const std::size_t jm = grid.idx(i, grid.wrap_y(j - 1), k);
                const double dpdy = (grid.p[jp] - grid.p[jm]) * inv_2dy;

                grid.du[c] += -inv_rho * dpdx;
                grid.dv[c] += -inv_rho * dpdy;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// coriolis: du += f*v,  dv += -f*u
// ---------------------------------------------------------------------------
void coriolis(Grid& grid, double latitude_rad) {
    const double f = atm::coriolis_parameter(latitude_rad);
    const std::size_t N = grid.total();

    for (std::size_t c = 0; c < N; ++c) {
        if (grid.solid[c]) continue;
        grid.du[c] +=  f * grid.v[c];
        grid.dv[c] += -f * grid.u[c];
    }
}

// ---------------------------------------------------------------------------
// Helper: Laplacian diffusion for a single scalar field.
// tendency += K * nabla^2(field)
// ---------------------------------------------------------------------------
static void diffuse_field(Grid& grid,
                          const std::vector<double>& field,
                          std::vector<double>& tendency,
                          double K) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;
    const double inv_dx2 = 1.0 / (grid.dx * grid.dx);
    const double inv_dy2 = 1.0 / (grid.dy * grid.dy);

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const std::size_t c = grid.idx(i, j, k);
                if (grid.solid[c]) continue;

                const double f_c = field[c];

                // Horizontal Laplacian: d2f/dx2 + d2f/dy2
                const double f_xp = field[grid.idx(grid.wrap_x(i + 1), j, k)];
                const double f_xm = field[grid.idx(grid.wrap_x(i - 1), j, k)];
                const double d2fdx2 = (f_xp - 2.0 * f_c + f_xm) * inv_dx2;

                const double f_yp = field[grid.idx(i, grid.wrap_y(j + 1), k)];
                const double f_ym = field[grid.idx(i, grid.wrap_y(j - 1), k)];
                const double d2fdy2 = (f_yp - 2.0 * f_c + f_ym) * inv_dy2;

                // Vertical Laplacian: d2f/dz2 (non-uniform spacing)
                double d2fdz2 = 0.0;
                const int kp = grid.clamp_z(k + 1);
                const int km = grid.clamp_z(k - 1);
                if (kp != k && km != k) {
                    // Interior: second derivative on non-uniform grid
                    const double dz_below = 0.5 * (grid.dz[km] + grid.dz[k]);
                    const double dz_above = 0.5 * (grid.dz[k] + grid.dz[kp]);
                    const double dz_avg = 0.5 * (dz_below + dz_above);
                    const double f_zp = field[grid.idx(i, j, kp)];
                    const double f_zm = field[grid.idx(i, j, km)];
                    d2fdz2 = ((f_zp - f_c) / dz_above - (f_c - f_zm) / dz_below) / dz_avg;
                } else if (kp != k) {
                    // Bottom boundary: zero-gradient (Neumann) condition
                    const double dz_above = 0.5 * (grid.dz[k] + grid.dz[kp]);
                    const double f_zp = field[grid.idx(i, j, kp)];
                    d2fdz2 = (f_zp - f_c) / (dz_above * dz_above);
                } else if (km != k) {
                    // Top boundary: zero-gradient (Neumann) condition
                    const double dz_below = 0.5 * (grid.dz[km] + grid.dz[k]);
                    const double f_zm = field[grid.idx(i, j, km)];
                    d2fdz2 = (f_zm - f_c) / (dz_below * dz_below);
                }
                // else single-level grid: d2fdz2 remains 0

                tendency[c] += K * (d2fdx2 + d2fdy2 + d2fdz2);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// diffusion: apply Laplacian diffusion to momentum (Km) and scalars (Kh)
// ---------------------------------------------------------------------------
void diffusion(Grid& grid, double Km, double Kh) {
    // Momentum diffusion
    diffuse_field(grid, grid.u, grid.du, Km);
    diffuse_field(grid, grid.v, grid.dv, Km);
    diffuse_field(grid, grid.w, grid.dw, Km);

    // Scalar (thermal + moisture) diffusion
    diffuse_field(grid, grid.theta, grid.dtheta, Kh);
    diffuse_field(grid, grid.qv,    grid.dqv,    Kh);
    diffuse_field(grid, grid.qc,    grid.dqc,    Kh);
    diffuse_field(grid, grid.qr,    grid.dqr,    Kh);
}

// ---------------------------------------------------------------------------
// apply_tendencies: forward Euler update for all prognostic variables.
// Moisture mixing ratios are clamped to non-negative values.
// ---------------------------------------------------------------------------
void apply_tendencies(Grid& grid, double dt) {
    const std::size_t N = grid.total();

    for (std::size_t c = 0; c < N; ++c) {
        if (grid.solid[c]) continue;

        grid.u[c]     += grid.du[c]     * dt;
        grid.v[c]     += grid.dv[c]     * dt;
        grid.w[c]     += grid.dw[c]     * dt;
        grid.theta[c] += grid.dtheta[c] * dt;
        grid.qv[c]    += grid.dqv[c]    * dt;
        grid.qc[c]    += grid.dqc[c]    * dt;
        grid.qr[c]    += grid.dqr[c]    * dt;

        // Moisture cannot be negative
        grid.qv[c] = std::max(grid.qv[c], 0.0);
        grid.qc[c] = std::max(grid.qc[c], 0.0);
        grid.qr[c] = std::max(grid.qr[c], 0.0);
    }
}

// ---------------------------------------------------------------------------
// sponge_layer: Rayleigh damping in the top sponge_depth levels.
// Relaxes prognostic variables toward the horizontally averaged state at each
// level. The damping coefficient increases linearly from zero at the base of
// the sponge to a maximum at the domain top.
// ---------------------------------------------------------------------------
void sponge_layer(Grid& grid, int sponge_depth) {
    if (sponge_depth <= 0) return;

    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;
    const int k_start = Nz - sponge_depth; // first sponge level
    if (k_start < 0) return;

    // Maximum inverse relaxation timescale [1/s]
    const double alpha_max = 0.05;

    for (int k = k_start; k < Nz; ++k) {
        // Compute horizontal mean at this level (reference state)
        double u_mean     = 0.0, v_mean     = 0.0, w_mean     = 0.0;
        double theta_mean = 0.0;
        double qv_mean    = 0.0, qc_mean    = 0.0, qr_mean    = 0.0;
        int count = 0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const std::size_t c = grid.idx(i, j, k);
                if (grid.solid[c]) continue;
                u_mean     += grid.u[c];
                v_mean     += grid.v[c];
                w_mean     += grid.w[c];
                theta_mean += grid.theta[c];
                qv_mean    += grid.qv[c];
                qc_mean    += grid.qc[c];
                qr_mean    += grid.qr[c];
                ++count;
            }
        }

        if (count == 0) continue;
        const double inv_count = 1.0 / static_cast<double>(count);
        u_mean     *= inv_count;
        v_mean     *= inv_count;
        w_mean     *= inv_count;
        theta_mean *= inv_count;
        qv_mean    *= inv_count;
        qc_mean    *= inv_count;
        qr_mean    *= inv_count;

        // Damping coefficient ramps linearly from 0 at k_start to alpha_max at Nz-1
        const double alpha = alpha_max
            * static_cast<double>(k - k_start + 1)
            / static_cast<double>(sponge_depth);

        // Apply Rayleigh relaxation to tendencies
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const std::size_t c = grid.idx(i, j, k);
                if (grid.solid[c]) continue;

                grid.du[c]     += -alpha * (grid.u[c]     - u_mean);
                grid.dv[c]     += -alpha * (grid.v[c]     - v_mean);
                grid.dw[c]     += -alpha * (grid.w[c]     - w_mean);
                grid.dtheta[c] += -alpha * (grid.theta[c] - theta_mean);
                grid.dqv[c]    += -alpha * (grid.qv[c]    - qv_mean);
                grid.dqc[c]    += -alpha * (grid.qc[c]    - qc_mean);
                grid.dqr[c]    += -alpha * (grid.qr[c]    - qr_mean);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// enforce_continuity: diagnose w from the anelastic continuity equation.
//   div(rho * V) = 0  =>  dw/dz = -(du/dx + dv/dy)  (Boussinesq)
// Integrates upward from w=0 at k=0 (surface no-penetration BC).
// Buoyancy-driven vertical acceleration is already present in the w field
// from apply_tendencies; this function blends with the continuity-diagnosed
// w to ensure mass conservation while preserving buoyancy forcing.
// ---------------------------------------------------------------------------
void enforce_continuity(Grid& grid) {
    const int Nx = grid.Nx;
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;
    const double inv_2dx = 1.0 / (2.0 * grid.dx);
    const double inv_2dy = 1.0 / (2.0 * grid.dy);

    // Blending factor: how much of the continuity-diagnosed w to use
    // 1.0 = fully diagnosed (pure hydrostatic), 0.0 = fully prognostic
    constexpr double blend = 0.8;

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double w_cont = 0.0; // w at surface = 0

            for (int k = 0; k < Nz; ++k) {
                std::size_t c = grid.idx(i, j, k);
                if (grid.solid[c]) {
                    w_cont = 0.0;
                    continue;
                }

                // Compute horizontal divergence at this cell
                std::size_t ip = grid.idx(grid.wrap_x(i + 1), j, k);
                std::size_t im = grid.idx(grid.wrap_x(i - 1), j, k);
                std::size_t jp = grid.idx(i, grid.wrap_y(j + 1), k);
                std::size_t jm = grid.idx(i, grid.wrap_y(j - 1), k);

                double dudx = (grid.u[ip] - grid.u[im]) * inv_2dx;
                double dvdy = (grid.v[jp] - grid.v[jm]) * inv_2dy;
                double div_h = dudx + dvdy;

                // Integrate continuity upward
                w_cont = w_cont - div_h * grid.dz[k];

                // Blend continuity-diagnosed w with prognostic w (from buoyancy)
                grid.w[c] = blend * w_cont + (1.0 - blend) * grid.w[c];
            }
        }
    }

    // Enforce w=0 at top boundary
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            std::size_t c = grid.idx(i, j, Nz - 1);
            grid.w[c] = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// clamp_velocities: limit wind speeds to physically reasonable bounds
// ---------------------------------------------------------------------------
void clamp_velocities(Grid& grid, double max_h, double max_v) {
    for (std::size_t c = 0; c < grid.total(); ++c) {
        if (grid.solid[c]) continue;
        grid.u[c] = std::max(-max_h, std::min(max_h, grid.u[c]));
        grid.v[c] = std::max(-max_h, std::min(max_h, grid.v[c]));
        grid.w[c] = std::max(-max_v, std::min(max_v, grid.w[c]));
    }
}

} // namespace physics
