#pragma once
#include <string>

struct Config {
    // Grid dimensions
    int nx = 100;
    int ny = 100;
    int nz = 60;

    // Grid spacing [m]
    double dx = 1000.0;
    double dy = 1000.0;
    double dz_base = 200.0;

    // Time
    double timestep = 2.0;
    int run_steps = 3600;
    int output_interval = 100;

    // Physics
    double latitude = -27.5;    // degrees, for Coriolis
    double Km = 500.0;          // momentum diffusion [m^2/s]
    double Kh = 500.0;          // scalar diffusion [m^2/s]
    int sponge_levels = 5;
    double solar_angle = 0.785; // solar zenith angle [rad]

    // Terrain
    std::string terrain_file = "";
    double terrain_max_height = 1500.0;

    // Data ingestion
    bool use_gfs = false;
    std::string gfs_file = "";
    std::string radiosonde_file = "";
    std::string aws_file = "";

    // Output
    std::string output_dir = "output/";
    bool output_vtk = true;
    bool output_csv = true;

    // Initial conditions
    double init_wind_u = 5.0;
    double init_wind_v = 0.0;
    double init_theta_surface = 300.0;
    double init_qv_surface = 0.012;

    // Warm bubble perturbation for convection tests
    bool add_warm_bubble = false;
    double bubble_dtheta = 3.0;
    double bubble_radius = 5000.0;
};

Config load_config(const std::string& path);
