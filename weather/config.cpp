#include "config.hpp"
#include <fstream>
#include <regex>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {
std::string trim(const std::string& s) {
    const auto start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    const auto end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}
}

Config load_config(const std::string& path) {
    Config cfg;
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Failed to open config: " + path);
    }

    // Map of floating point fields for easy assignment
    std::unordered_map<std::string, double*> float_fields{
        {"timestep", &cfg.timestep}
    };

    std::string line;
    std::regex domain_re(R"(domain_size\s*=\s*\[(\d+),\s*(\d+),\s*(\d+)\])");
    std::regex float_re(R"((\w+)\s*=\s*([0-9]*\.?[0-9]+))");
    std::regex bool_re(R"((\w+)\s*=\s*(true|false))");
    std::regex string_re(R"regex((\w+)\s*=\s*"([^"]*)")regex");
    std::regex int_re(R"((\w+)\s*=\s*(\d+))");
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::smatch m;
        if (std::regex_match(line, m, domain_re)) {
            cfg.domain_size = {std::stoi(m[1]), std::stoi(m[2]), std::stoi(m[3])};
        } else if (std::regex_match(line, m, int_re) && m[1] == "run_steps") {
            cfg.run_steps = std::stoi(m[2]);
        } else if (std::regex_match(line, m, float_re)) {
            auto it = float_fields.find(m[1]);
            if (it != float_fields.end()) {
                *(it->second) = std::stod(m[2]);
            }
        } else if (std::regex_match(line, m, bool_re) && m[1] == "use_gfs") {
            cfg.use_gfs = (m[2] == "true");
        } else if (std::regex_match(line, m, string_re)) {
            std::string key = m[1];
            std::string val = m[2];
            if (key == "gfs_file") cfg.gfs_file = val;
            else if (key == "radiosonde_file") cfg.radiosonde_file = val;
            else if (key == "aws_file") cfg.aws_file = val;
            else if (key == "terrain_file") cfg.terrain_file = val;
            else if (key == "output_dir") cfg.output_dir = val;
        }
    }
    return cfg;
}
