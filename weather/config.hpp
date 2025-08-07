#pragma once
#include <array>
#include <string>

struct Config {
    std::array<int,3> domain_size{};
    double timestep{};
    bool use_gfs{};
    std::string gfs_file;
    std::string radiosonde_file;
    std::string aws_file;
    std::string terrain_file;
    std::string output_dir;
    int run_steps{};
};

Config load_config(const std::string& path);
