#include "radiosonde.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace io {

std::vector<RadiosondeRecord> ingest_radiosonde(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Failed to open radiosonde file: " + path);
    }
    std::vector<RadiosondeRecord> data;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        RadiosondeRecord rec{};
        if (iss >> rec.pressure >> rec.temperature >> rec.humidity) {
            data.push_back(rec);
        }
    }
    return data;
}

} // namespace io

