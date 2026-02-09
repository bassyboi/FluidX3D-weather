#include "grid.hpp"
#include "constants.hpp"
#include <cmath>
#include <algorithm>

Grid::Grid(int nx, int ny, int nz, double dx_m, double dy_m, double dz_base)
    : Nx(nx), Ny(ny), Nz(nz), dx(dx_m), dy(dy_m) {
    init_vertical_levels(dz_base);
    allocate();
}

void Grid::init_vertical_levels(double dz_base) {
    dz.resize(Nz);
    z.resize(Nz);
    double z_acc = 0.0;
    for (int k = 0; k < Nz; ++k) {
        double stretch = 1.0 + 0.02 * k; // gentle stretching: finer near surface
        dz[k] = dz_base * stretch;
        z[k] = z_acc + dz[k] * 0.5; // cell center height
        z_acc += dz[k];
    }
}

void Grid::allocate() {
    std::size_t N = total();
    std::size_t N2 = (std::size_t)Nx * Ny;

    u.assign(N, 0.0);
    v.assign(N, 0.0);
    w.assign(N, 0.0);
    theta.assign(N, 300.0);
    qv.assign(N, 0.005);
    qc.assign(N, 0.0);
    qr.assign(N, 0.0);

    p.assign(N, atm::p0);
    rho.assign(N, atm::rho0);
    T.assign(N, atm::T0);

    terrain_height.assign(N2, 0.0);
    solid.assign(N, false);

    du.assign(N, 0.0);
    dv.assign(N, 0.0);
    dw.assign(N, 0.0);
    dtheta.assign(N, 0.0);
    dqv.assign(N, 0.0);
    dqc.assign(N, 0.0);
    dqr.assign(N, 0.0);
}

void Grid::compute_diagnostics() {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            // Set pressure at top level using barometric formula
            int k_top = Nz - 1;
            std::size_t idx_top = idx(i, j, k_top);
            double p_top = atm::p0 * std::pow(
                std::max(0.0, 1.0 - atm::g * z[k_top] / (atm::cp_d * atm::T0)),
                atm::cp_d / atm::Rd);
            p[idx_top] = std::max(p_top, 5000.0);

            T[idx_top] = atm::temperature_from_theta(theta[idx_top], p[idx_top]);
            double Tv = atm::virtual_temperature(T[idx_top], std::max(qv[idx_top], 0.0));
            rho[idx_top] = p[idx_top] / (atm::Rd * Tv);

            // Integrate downward hydrostatically
            for (int k = Nz - 2; k >= 0; --k) {
                std::size_t id = idx(i, j, k);
                std::size_t id_above = idx(i, j, k + 1);

                double dz_between = 0.5 * (dz[k] + dz[k + 1]);
                double theta_avg = 0.5 * (theta[id] + theta[id_above]);
                double qv_avg = 0.5 * (std::max(qv[id], 0.0) + std::max(qv[id_above], 0.0));

                double T_guess = atm::temperature_from_theta(theta_avg, p[id_above]);
                double Tv_guess = atm::virtual_temperature(T_guess, qv_avg);
                double rho_avg = p[id_above] / (atm::Rd * Tv_guess);

                p[id] = p[id_above] + rho_avg * atm::g * dz_between;
                T[id] = atm::temperature_from_theta(theta[id], p[id]);
                Tv = atm::virtual_temperature(T[id], std::max(qv[id], 0.0));
                rho[id] = p[id] / (atm::Rd * Tv);
            }
        }
    }
}

void Grid::zero_tendencies() {
    std::fill(du.begin(), du.end(), 0.0);
    std::fill(dv.begin(), dv.end(), 0.0);
    std::fill(dw.begin(), dw.end(), 0.0);
    std::fill(dtheta.begin(), dtheta.end(), 0.0);
    std::fill(dqv.begin(), dqv.end(), 0.0);
    std::fill(dqc.begin(), dqc.end(), 0.0);
    std::fill(dqr.begin(), dqr.end(), 0.0);
}
