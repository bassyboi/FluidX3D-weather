#include "config.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace {

std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::string unquote(const std::string& s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        return s.substr(1, s.size() - 2);
    return s;
}

} // anonymous namespace

Config load_config(const std::string& path) {
    Config cfg;
    std::ifstream file(path);
    if (!file) {
        std::cerr << "Warning: Could not open config file '" << path
                  << "', using defaults.\n";
        return cfg;
    }

    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));

        // Handle array syntax: domain_size = [100, 100, 60]
        if (key == "domain_size") {
            std::string nums = val;
            nums.erase(std::remove(nums.begin(), nums.end(), '['), nums.end());
            nums.erase(std::remove(nums.begin(), nums.end(), ']'), nums.end());
            std::istringstream iss(nums);
            char comma;
            iss >> cfg.nx >> comma >> cfg.ny >> comma >> cfg.nz;
            continue;
        }

        std::string sv = unquote(val);

        if      (key == "nx")                cfg.nx = std::stoi(val);
        else if (key == "ny")                cfg.ny = std::stoi(val);
        else if (key == "nz")                cfg.nz = std::stoi(val);
        else if (key == "dx")                cfg.dx = std::stod(val);
        else if (key == "dy")                cfg.dy = std::stod(val);
        else if (key == "dz_base")           cfg.dz_base = std::stod(val);
        else if (key == "timestep")          cfg.timestep = std::stod(val);
        else if (key == "run_steps")         cfg.run_steps = std::stoi(val);
        else if (key == "output_interval")   cfg.output_interval = std::stoi(val);
        else if (key == "latitude")          cfg.latitude = std::stod(val);
        else if (key == "Km")               cfg.Km = std::stod(val);
        else if (key == "Kh")               cfg.Kh = std::stod(val);
        else if (key == "sponge_levels")     cfg.sponge_levels = std::stoi(val);
        else if (key == "solar_angle")       cfg.solar_angle = std::stod(val);
        else if (key == "terrain_file")      cfg.terrain_file = sv;
        else if (key == "terrain_max_height") cfg.terrain_max_height = std::stod(val);
        else if (key == "use_gfs")           cfg.use_gfs = (val == "true" || val == "1");
        else if (key == "gfs_file")          cfg.gfs_file = sv;
        else if (key == "radiosonde_file")   cfg.radiosonde_file = sv;
        else if (key == "aws_file")          cfg.aws_file = sv;
        else if (key == "output_dir")        cfg.output_dir = sv;
        else if (key == "output_vtk")        cfg.output_vtk = (val == "true" || val == "1");
        else if (key == "output_csv")        cfg.output_csv = (val == "true" || val == "1");
        else if (key == "init_wind_u")       cfg.init_wind_u = std::stod(val);
        else if (key == "init_wind_v")       cfg.init_wind_v = std::stod(val);
        else if (key == "init_theta_surface") cfg.init_theta_surface = std::stod(val);
        else if (key == "init_qv_surface")   cfg.init_qv_surface = std::stod(val);
        else if (key == "add_warm_bubble")   cfg.add_warm_bubble = (val == "true" || val == "1");
        else if (key == "bubble_dtheta")     cfg.bubble_dtheta = std::stod(val);
        else if (key == "bubble_radius")     cfg.bubble_radius = std::stod(val);
        else {
            std::cerr << "Warning: Unknown config key '" << key << "'\n";
        }
    }

    return cfg;
}
