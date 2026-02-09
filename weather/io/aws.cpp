#include "aws.hpp"
#include "../constants.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

namespace io {

// A single automatic weather station record
struct AWSRecord {
    double time_s;       // time [s]
    double temperature_C;// surface temperature [C]
    double humidity_pct; // relative humidity [%]
    double pressure_hPa; // surface pressure [hPa]
    double wind_speed_ms;// wind speed [m/s]
    double wind_dir_deg; // wind direction [degrees from north]
};

void ingest_aws(Grid& grid, const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[aws] Warning: could not open '" << path << "'\n";
        return;
    }

    // Read all records, but we only use the first one for initial conditions
    AWSRecord first_rec;
    bool found = false;
    std::string line;

    // Skip header line
    if (std::getline(file, line)) {
        // Check if it looks like data
        if (!line.empty() && (std::isdigit(static_cast<unsigned char>(line[0])) ||
                              line[0] == '-' || line[0] == '+')) {
            std::istringstream iss(line);
            char comma;
            if (iss >> first_rec.time_s >> comma >> first_rec.temperature_C >> comma
                    >> first_rec.humidity_pct >> comma >> first_rec.pressure_hPa >> comma
                    >> first_rec.wind_speed_ms >> comma >> first_rec.wind_dir_deg) {
                found = true;
            }
        }
    }

    if (!found) {
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            char comma;
            if (iss >> first_rec.time_s >> comma >> first_rec.temperature_C >> comma
                    >> first_rec.humidity_pct >> comma >> first_rec.pressure_hPa >> comma
                    >> first_rec.wind_speed_ms >> comma >> first_rec.wind_dir_deg) {
                found = true;
                break;
            }
        }
    }

    if (!found) {
        std::cerr << "[aws] No valid records found in '" << path << "'\n";
        return;
    }

    std::cout << "[aws] Using surface observation: T=" << first_rec.temperature_C
              << " C, RH=" << first_rec.humidity_pct
              << " %, P=" << first_rec.pressure_hPa
              << " hPa, wind=" << first_rec.wind_speed_ms << " m/s from "
              << first_rec.wind_dir_deg << " deg\n";

    // Convert to model variables
    double T_K    = first_rec.temperature_C + 273.15;
    double p_Pa   = first_rec.pressure_hPa * 100.0;
    double theta  = atm::potential_temperature(T_K, p_Pa);

    // Convert relative humidity to mixing ratio
    double qvs = atm::saturation_mixing_ratio(p_Pa, T_K);
    double qv  = (first_rec.humidity_pct / 100.0) * qvs;
    qv = std::max(qv, 0.0);

    // Convert wind direction and speed to u, v components
    double dir_rad = first_rec.wind_dir_deg * M_PI / 180.0;
    double u_wind = -first_rec.wind_speed_ms * std::sin(dir_rad);
    double v_wind = -first_rec.wind_speed_ms * std::cos(dir_rad);

    double rho_val = atm::density(p_Pa, T_K);

    // Apply surface conditions at k=0 across the entire domain
    int k = 0;
    for (int j = 0; j < grid.Ny; ++j) {
        for (int i = 0; i < grid.Nx; ++i) {
            std::size_t id = grid.idx(i, j, k);
            if (grid.solid[id]) continue;

            grid.theta[id] = theta;
            grid.qv[id]    = qv;
            grid.u[id]     = u_wind;
            grid.v[id]     = v_wind;
            grid.w[id]     = 0.0;
            grid.p[id]     = p_Pa;
            grid.T[id]     = T_K;
            grid.rho[id]   = rho_val;
        }
    }

    std::cout << "[aws] Applied surface boundary conditions at k=0.\n";
}

} // namespace io
