#include "config.hpp"
#include "grid.hpp"
#include "constants.hpp"
#include "terrain.hpp"
#include "physics/solver.hpp"
#include "physics/thermal.hpp"
#include "physics/cloud.hpp"
#include "io/gfs.hpp"
#include "io/radiosonde.hpp"
#include "io/aws.hpp"
#include "io/output.hpp"
#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>

// Initialize atmospheric state with standard atmosphere + optional perturbation
static void initialize_atmosphere(Grid& grid, const Config& cfg) {
    double theta_sfc = cfg.init_theta_surface;
    double qv_sfc = cfg.init_qv_surface;
    double lapse = 6.5e-3; // K/m tropospheric lapse rate

    for (int k = 0; k < grid.Nz; ++k) {
        double z_m = grid.z[k];
        // Standard atmosphere temperature profile
        double T_std = atm::T0 - lapse * z_m;
        if (T_std < 200.0) T_std = 200.0; // tropopause minimum

        // Pressure from barometric formula
        double p_std = atm::p0 * std::pow(
            std::max(0.001, 1.0 - lapse * z_m / atm::T0),
            atm::g / (atm::Rd * lapse));

        // Potential temperature
        double theta_std = T_std * std::pow(atm::p0 / p_std, atm::kappa);

        // Scale surface values to this level
        // Add slight stable stratification above surface
        double theta_val = theta_sfc + 3.0e-3 * z_m; // ~3 K/km increase in theta
        if (theta_val < theta_std) theta_val = theta_std;

        // Moisture decreases exponentially with height
        double qv_val = qv_sfc * std::exp(-z_m / 2500.0);
        qv_val = std::max(qv_val, 1.0e-6);

        for (int i = 0; i < grid.Nx; ++i) {
            for (int j = 0; j < grid.Ny; ++j) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;

                grid.theta[id] = theta_val;
                grid.qv[id] = qv_val;
                grid.u[id] = cfg.init_wind_u;
                grid.v[id] = cfg.init_wind_v;
                grid.w[id] = 0.0;
                grid.qc[id] = 0.0;
                grid.qr[id] = 0.0;
            }
        }
    }

    // Add warm bubble perturbation for convection testing
    if (cfg.add_warm_bubble) {
        double cx = grid.Nx * grid.dx * 0.5;
        double cy = grid.Ny * grid.dy * 0.5;
        double cz = 1500.0; // bubble center at 1.5 km
        double R = cfg.bubble_radius;

        for (int k = 0; k < grid.Nz; ++k) {
            for (int j = 0; j < grid.Ny; ++j) {
                for (int i = 0; i < grid.Nx; ++i) {
                    std::size_t id = grid.idx(i, j, k);
                    if (grid.solid[id]) continue;

                    double x = (i + 0.5) * grid.dx;
                    double y = (j + 0.5) * grid.dy;
                    double z = grid.z[k];
                    double r = std::sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));

                    if (r < R) {
                        double factor = std::cos(0.5 * M_PI * r / R);
                        grid.theta[id] += cfg.bubble_dtheta * factor * factor;
                    }
                }
            }
        }
        std::cout << "  Added warm bubble: dtheta=" << cfg.bubble_dtheta
                  << " K, radius=" << R/1000.0 << " km, center=("
                  << cx/1000.0 << "," << cy/1000.0 << "," << cz/1000.0 << ") km\n";
    }
}

// Print simulation progress
static void print_progress(int step, int total, const Grid& grid, double elapsed_s) {
    // Find max wind speed and temperature range
    double max_wind = 0.0;
    double min_T = 1e10, max_T = -1e10;
    double total_qc = 0.0, total_qr = 0.0;

    for (std::size_t n = 0; n < grid.total(); ++n) {
        if (grid.solid[n]) continue;
        double wind = std::sqrt(grid.u[n]*grid.u[n] + grid.v[n]*grid.v[n] + grid.w[n]*grid.w[n]);
        if (wind > max_wind) max_wind = wind;
        if (grid.T[n] < min_T) min_T = grid.T[n];
        if (grid.T[n] > max_T) max_T = grid.T[n];
        total_qc += grid.qc[n];
        total_qr += grid.qr[n];
    }

    double pct = 100.0 * step / total;
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Step " << step << "/" << total << " (" << pct << "%)";
    std::cout << "  |V|max=" << std::setprecision(2) << max_wind << " m/s";
    std::cout << "  T=[" << std::setprecision(1) << (min_T-273.15)
              << "," << (max_T-273.15) << "] C";
    if (total_qc > 0.0 || total_qr > 0.0)
        std::cout << "  cloud=" << std::scientific << std::setprecision(2) << total_qc
                  << "  rain=" << total_qr;
    std::cout << "  elapsed=" << std::fixed << std::setprecision(1) << elapsed_s << "s\n";
}

int main(int argc, char* argv[]) {
    std::cout << "========================================\n";
    std::cout << " FluidX3D-Weather Prediction Model\n";
    std::cout << "========================================\n\n";

    // Load configuration
    std::string config_path = "config.toml";
    if (argc > 1) config_path = argv[1];
    std::cout << "Loading config from: " << config_path << "\n";
    Config cfg = load_config(config_path);

    // Create grid
    std::cout << "Creating grid: " << cfg.nx << "x" << cfg.ny << "x" << cfg.nz
              << " (dx=" << cfg.dx << "m, dy=" << cfg.dy << "m, dz_base=" << cfg.dz_base << "m)\n";
    Grid grid(cfg.nx, cfg.ny, cfg.nz, cfg.dx, cfg.dy, cfg.dz_base);

    double domain_x_km = cfg.nx * cfg.dx / 1000.0;
    double domain_y_km = cfg.ny * cfg.dy / 1000.0;
    double domain_z_km = grid.z[grid.Nz-1] / 1000.0;
    std::cout << "Domain size: " << domain_x_km << " x " << domain_y_km
              << " x " << std::fixed << std::setprecision(1) << domain_z_km << " km\n";

    // Load/generate terrain
    std::cout << "Setting up terrain...\n";
    if (!cfg.terrain_file.empty()) {
        load_terrain(grid, cfg.terrain_file);
    } else {
        generate_terrain(grid, cfg.terrain_max_height);
    }

    // Initialize atmosphere
    std::cout << "Initializing atmosphere...\n";
    initialize_atmosphere(grid, cfg);

    // Ingest observational data if available
    if (cfg.use_gfs && !cfg.gfs_file.empty()) {
        std::cout << "Loading GFS data: " << cfg.gfs_file << "\n";
        io::load_gfs(grid, cfg.gfs_file);
    }
    if (!cfg.radiosonde_file.empty()) {
        std::cout << "Loading radiosonde data: " << cfg.radiosonde_file << "\n";
        io::ingest_radiosonde(grid, cfg.radiosonde_file);
    }
    if (!cfg.aws_file.empty()) {
        std::cout << "Loading AWS data: " << cfg.aws_file << "\n";
        io::ingest_aws(grid, cfg.aws_file);
    }

    // Compute initial diagnostics
    grid.compute_diagnostics();

    // Set up output
    io::ensure_directory(cfg.output_dir);

    // Convert latitude to radians for Coriolis
    double lat_rad = cfg.latitude * M_PI / 180.0;

    std::cout << "\nSimulation parameters:\n";
    std::cout << "  dt = " << cfg.timestep << " s\n";
    std::cout << "  steps = " << cfg.run_steps << "\n";
    std::cout << "  total time = " << cfg.run_steps * cfg.timestep / 3600.0 << " hours\n";
    std::cout << "  latitude = " << cfg.latitude << " deg (f = "
              << std::scientific << std::setprecision(3) << atm::coriolis_parameter(lat_rad) << " /s)\n";
    std::cout << "  Km = " << std::fixed << cfg.Km << " m^2/s, Kh = " << cfg.Kh << " m^2/s\n";
    std::cout << "\nStarting simulation...\n";

    auto t_start = std::chrono::steady_clock::now();

    // Write initial state
    if (cfg.output_vtk) io::write_vtk(grid, cfg.output_dir, 0);
    if (cfg.output_csv) io::write_csv_summary(grid, cfg.output_dir, 0);

    // Main time loop
    for (int step = 1; step <= cfg.run_steps; ++step) {
        double dt = cfg.timestep;

        // 1. Zero tendencies
        grid.zero_tendencies();

        // 2. Dynamics: advection
        physics::advect(grid, dt);

        // 3. Pressure gradient force
        physics::pressure_gradient(grid);

        // 4. Coriolis force
        physics::coriolis(grid, lat_rad);

        // 5. Thermal buoyancy
        physics::apply_thermal_buoyancy(grid);

        // 6. Radiation (Newtonian cooling)
        physics::apply_radiation(grid, dt);

        // 7. Surface heat flux
        physics::surface_heat_flux(grid, cfg.solar_angle, dt);

        // 8. Diffusion (momentum + scalar)
        physics::diffusion(grid, cfg.Km, cfg.Kh);

        // 9. Sponge layer damping at domain top
        physics::sponge_layer(grid, cfg.sponge_levels);

        // 10. Apply all accumulated tendencies (forward Euler)
        physics::apply_tendencies(grid, dt);

        // 11. Enforce mass continuity (diagnose w from horizontal divergence)
        physics::enforce_continuity(grid);

        // 12. Clamp velocities for numerical safety
        physics::clamp_velocities(grid, 80.0, 40.0);

        // 13. Cloud microphysics (applied directly, not via tendencies)
        physics::update_cloud_microphysics(grid, dt);

        // 14. Recompute diagnostic variables (p, T, rho)
        grid.compute_diagnostics();

        // Progress and output
        if (step % cfg.output_interval == 0 || step == cfg.run_steps) {
            auto t_now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(t_now - t_start).count();
            print_progress(step, cfg.run_steps, grid, elapsed);

            if (cfg.output_vtk) io::write_vtk(grid, cfg.output_dir, step);
            if (cfg.output_csv) io::write_csv_summary(grid, cfg.output_dir, step);
        }
    }

    auto t_end = std::chrono::steady_clock::now();
    double total_time = std::chrono::duration<double>(t_end - t_start).count();

    std::cout << "\n========================================\n";
    std::cout << " Simulation complete!\n";
    std::cout << " Total wall time: " << std::fixed << std::setprecision(1) << total_time << " s\n";
    std::cout << " Simulated time: " << cfg.run_steps * cfg.timestep / 3600.0 << " hours\n";
    std::cout << " Output in: " << cfg.output_dir << "\n";
    std::cout << "========================================\n";

    return 0;
}
