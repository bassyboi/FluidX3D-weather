#pragma once
#include <vector>
#include <cstddef>

// 3D atmospheric grid with all prognostic and diagnostic variables
struct Grid {
    int Nx, Ny, Nz;        // grid dimensions
    double dx, dy;          // horizontal grid spacing [m]
    std::vector<double> dz; // vertical grid spacing per level [m]
    std::vector<double> z;  // height of each level [m]

    // Prognostic variables (evolved in time)
    std::vector<double> u;      // x-wind component [m/s]
    std::vector<double> v;      // y-wind component [m/s]
    std::vector<double> w;      // z-wind component [m/s]
    std::vector<double> theta;  // potential temperature [K]
    std::vector<double> qv;     // water vapor mixing ratio [kg/kg]
    std::vector<double> qc;     // cloud water mixing ratio [kg/kg]
    std::vector<double> qr;     // rain water mixing ratio [kg/kg]

    // Diagnostic variables (derived each step)
    std::vector<double> p;      // pressure [Pa]
    std::vector<double> rho;    // density [kg/m^3]
    std::vector<double> T;      // temperature [K]

    // Terrain
    std::vector<double> terrain_height; // terrain elevation [m], Nx*Ny
    std::vector<bool>   solid;          // true if cell is blocked by terrain

    // Tendencies (time derivatives for each prognostic variable)
    std::vector<double> du, dv, dw, dtheta, dqv, dqc, dqr;

    Grid() : Nx(0), Ny(0), Nz(0), dx(0), dy(0) {}
    Grid(int nx, int ny, int nz, double dx_m, double dy_m, double dz_base);

    void allocate();
    void init_vertical_levels(double dz_base);
    void compute_diagnostics(); // compute p, rho, T from theta, qv

    std::size_t total() const { return (std::size_t)Nx * Ny * Nz; }
    std::size_t idx(int x, int y, int zz) const {
        return (std::size_t)zz * Nx * Ny + (std::size_t)y * Nx + x;
    }

    int wrap_x(int x) const { return ((x % Nx) + Nx) % Nx; }
    int wrap_y(int y) const { return ((y % Ny) + Ny) % Ny; }
    int clamp_z(int zz) const { return zz < 0 ? 0 : (zz >= Nz ? Nz - 1 : zz); }

    void zero_tendencies();
};
