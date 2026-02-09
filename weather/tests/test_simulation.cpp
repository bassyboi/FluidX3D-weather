#include "../config.hpp"
#include "../grid.hpp"
#include "../constants.hpp"
#include "../terrain.hpp"
#include "../physics/solver.hpp"
#include "../physics/thermal.hpp"
#include "../physics/cloud.hpp"
#include "../io/gfs.hpp"
#include "../io/radiosonde.hpp"
#include "../io/aws.hpp"
#include "../io/output.hpp"
#include <cassert>
#include <iostream>
#include <cmath>
#include <fstream>

static void test_constants() {
    // Test saturation vapor pressure at 0C
    double es0 = atm::saturation_vapor_pressure(0.0);
    assert(std::abs(es0 - 611.2) < 1.0);

    // Test at 20C (~2338 Pa)
    double es20 = atm::saturation_vapor_pressure(20.0);
    assert(es20 > 2000.0 && es20 < 2800.0);

    // Test potential temperature
    double theta = atm::potential_temperature(288.15, 100000.0);
    assert(std::abs(theta - 288.15) < 0.1); // at p0, theta == T

    double theta2 = atm::potential_temperature(250.0, 50000.0);
    assert(theta2 > 250.0); // theta > T when p < p0

    // Test Coriolis parameter
    double f = atm::coriolis_parameter(M_PI / 4.0); // 45 degrees
    assert(std::abs(f - 1.0313e-4) < 1e-6);

    std::cout << "  constants: OK\n";
}

static void test_grid() {
    Grid grid(10, 10, 5, 1000.0, 1000.0, 200.0);

    assert(grid.Nx == 10);
    assert(grid.Ny == 10);
    assert(grid.Nz == 5);
    assert(grid.total() == 500);

    // Check vertical levels
    assert(grid.z[0] > 0.0);
    assert(grid.z[0] < grid.z[1]);
    assert(grid.dz[0] > 0.0);

    // Check index function
    assert(grid.idx(0, 0, 0) == 0);
    assert(grid.idx(1, 0, 0) == 1);
    assert(grid.idx(0, 1, 0) == 10);
    assert(grid.idx(0, 0, 1) == 100);

    // Check wrapping
    assert(grid.wrap_x(-1) == 9);
    assert(grid.wrap_x(10) == 0);
    assert(grid.wrap_y(-1) == 9);
    assert(grid.clamp_z(-1) == 0);
    assert(grid.clamp_z(5) == 4);

    // Check diagnostics computation
    grid.compute_diagnostics();
    for (int k = 0; k < grid.Nz; ++k) {
        std::size_t id = grid.idx(5, 5, k);
        assert(grid.p[id] > 0.0);
        assert(grid.T[id] > 100.0 && grid.T[id] < 400.0);
        assert(grid.rho[id] > 0.0 && grid.rho[id] < 3.0);
    }

    // Pressure should decrease with height
    assert(grid.p[grid.idx(5, 5, 0)] > grid.p[grid.idx(5, 5, 4)]);

    std::cout << "  grid: OK\n";
}

static void test_terrain() {
    Grid grid(20, 20, 10, 1000.0, 1000.0, 200.0);
    generate_terrain(grid, 500.0);

    // Some cells should be solid near terrain
    bool has_solid = false;
    bool has_fluid = false;
    for (std::size_t n = 0; n < grid.total(); ++n) {
        if (grid.solid[n]) has_solid = true;
        else has_fluid = true;
    }
    assert(has_solid); // terrain should create solid cells
    assert(has_fluid); // not everything should be solid

    // Solid cells should have zero velocity
    for (std::size_t n = 0; n < grid.total(); ++n) {
        if (grid.solid[n]) {
            assert(grid.u[n] == 0.0);
            assert(grid.v[n] == 0.0);
            assert(grid.w[n] == 0.0);
        }
    }

    std::cout << "  terrain: OK\n";
}

static void test_solver() {
    Grid grid(10, 10, 5, 1000.0, 1000.0, 200.0);

    // Set up a simple wind field
    for (std::size_t n = 0; n < grid.total(); ++n) {
        grid.u[n] = 5.0;
        grid.v[n] = 2.0;
        grid.w[n] = 0.0;
        grid.theta[n] = 300.0;
        grid.qv[n] = 0.01;
    }
    grid.compute_diagnostics();

    // Test advection
    grid.zero_tendencies();
    physics::advect(grid, 1.0);

    // Uniform field should have zero advection tendency
    for (std::size_t n = 0; n < grid.total(); ++n) {
        assert(std::abs(grid.dtheta[n]) < 1e-10);
    }

    // Test Coriolis
    grid.zero_tendencies();
    physics::coriolis(grid, M_PI / 6.0); // 30 degrees

    double f = atm::coriolis_parameter(M_PI / 6.0);
    std::size_t test_id = grid.idx(5, 5, 2);
    assert(std::abs(grid.du[test_id] - f * grid.v[test_id]) < 1e-10);
    assert(std::abs(grid.dv[test_id] - (-f * grid.u[test_id])) < 1e-10);

    // Test apply_tendencies
    double u_before = grid.u[test_id];
    grid.zero_tendencies();
    grid.du[test_id] = 1.0; // 1 m/s^2 acceleration
    physics::apply_tendencies(grid, 2.0);
    assert(std::abs(grid.u[test_id] - (u_before + 2.0)) < 1e-10);

    std::cout << "  solver: OK\n";
}

static void test_thermal() {
    Grid grid(10, 10, 5, 1000.0, 1000.0, 200.0);

    for (std::size_t n = 0; n < grid.total(); ++n) {
        grid.theta[n] = 300.0;
        grid.qv[n] = 0.01;
        grid.qc[n] = 0.0;
        grid.qr[n] = 0.0;
    }

    // Add warm perturbation at center
    std::size_t center = grid.idx(5, 5, 2);
    grid.theta[center] = 303.0;

    grid.compute_diagnostics();
    grid.zero_tendencies();
    physics::apply_thermal_buoyancy(grid);

    // Warm cell should have positive buoyancy (upward acceleration)
    assert(grid.dw[center] > 0.0);

    std::cout << "  thermal: OK\n";
}

static void test_cloud_microphysics() {
    Grid grid(5, 5, 3, 1000.0, 1000.0, 200.0);

    // Set up a supersaturated parcel
    for (std::size_t n = 0; n < grid.total(); ++n) {
        grid.theta[n] = 300.0;
        grid.qv[n] = 0.005;
        grid.qc[n] = 0.0;
        grid.qr[n] = 0.0;
    }
    grid.compute_diagnostics();

    // Make one cell highly supersaturated
    std::size_t test_id = grid.idx(2, 2, 1);
    double qvs = atm::saturation_mixing_ratio(grid.p[test_id], grid.T[test_id]);
    grid.qv[test_id] = qvs * 1.5; // 50% supersaturated

    double qv_before = grid.qv[test_id];
    double qc_before = grid.qc[test_id];

    physics::update_cloud_microphysics(grid, 1.0);

    // Should have condensed some water
    assert(grid.qv[test_id] < qv_before);
    assert(grid.qc[test_id] > qc_before);

    // Total water should be approximately conserved
    double total_before = qv_before + qc_before;
    double total_after = grid.qv[test_id] + grid.qc[test_id] + grid.qr[test_id];
    // Saturation adjustment changes theta which affects qvs,
    // so total water is approximately (not exactly) conserved
    assert(std::abs(total_before - total_after) / total_before < 0.1);

    std::cout << "  cloud microphysics: OK\n";
}

static void test_output() {
    Grid grid(5, 5, 3, 1000.0, 1000.0, 200.0);
    grid.compute_diagnostics();

    io::ensure_directory("test_output/");
    io::write_vtk(grid, "test_output/", 0);
    io::write_csv_summary(grid, "test_output/", 0);

    // Check files were created (by trying to open them)
    std::ifstream vtk("test_output/weather_000000.vtk");
    assert(vtk.good());
    std::ifstream csv("test_output/summary_000000.csv");
    assert(csv.good());

    std::cout << "  output: OK\n";
}

static void test_short_simulation() {
    Grid grid(10, 10, 5, 2000.0, 2000.0, 250.0);

    // Initialize
    for (std::size_t n = 0; n < grid.total(); ++n) {
        grid.theta[n] = 300.0 + 3.0e-3 * grid.z[n % grid.Nz];
        grid.qv[n] = 0.01 * std::exp(-grid.z[n % grid.Nz] / 2500.0);
        grid.u[n] = 5.0;
        grid.v[n] = 0.0;
        grid.w[n] = 0.0;
    }
    grid.compute_diagnostics();

    double dt = 2.0;
    int steps = 10;

    for (int step = 0; step < steps; ++step) {
        grid.zero_tendencies();
        physics::advect(grid, dt);
        physics::pressure_gradient(grid);
        physics::coriolis(grid, -0.48); // ~27.5S
        physics::apply_thermal_buoyancy(grid);
        physics::apply_radiation(grid, dt);
        physics::diffusion(grid, 500.0, 500.0);
        physics::sponge_layer(grid, 2);
        physics::apply_tendencies(grid, dt);
        physics::update_cloud_microphysics(grid, dt);
        grid.compute_diagnostics();
    }

    // Check simulation didn't blow up
    for (std::size_t n = 0; n < grid.total(); ++n) {
        assert(std::isfinite(grid.u[n]));
        assert(std::isfinite(grid.v[n]));
        assert(std::isfinite(grid.w[n]));
        assert(std::isfinite(grid.theta[n]));
        assert(std::isfinite(grid.T[n]));
        assert(std::isfinite(grid.p[n]));
        assert(grid.p[n] > 0.0);
        assert(grid.T[n] > 100.0 && grid.T[n] < 500.0);
        assert(grid.qv[n] >= 0.0);
        assert(grid.qc[n] >= 0.0);
        assert(grid.qr[n] >= 0.0);
    }

    std::cout << "  short simulation (10 steps): OK\n";
}

int main() {
    std::cout << "Running weather model tests...\n";

    test_constants();
    test_grid();
    test_terrain();
    test_solver();
    test_thermal();
    test_cloud_microphysics();
    test_output();
    test_short_simulation();

    std::cout << "\nAll tests passed!\n";
    return 0;
}
