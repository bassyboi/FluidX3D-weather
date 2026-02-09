#include "radiosonde.hpp"
#include "../constants.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace io {

// A single radiosonde observation record
struct SondeRecord {
    double p_hPa;       // pressure [hPa]
    double height_m;    // height [m]
    double T_C;         // temperature [C]
    double Td_C;        // dewpoint [C]
    double wind_dir;    // wind direction [degrees from north]
    double wind_speed;  // wind speed [knots]
};

// Convert knots to m/s
static constexpr double KNOTS_TO_MS = 0.5144;

// Convert wind direction and speed to u, v components
// Meteorological convention: wind FROM the given direction
static void wind_to_uv(double dir_deg, double speed_knots,
                        double& u_out, double& v_out) {
    double speed_ms = speed_knots * KNOTS_TO_MS;
    double dir_rad  = dir_deg * M_PI / 180.0;
    u_out = -speed_ms * std::sin(dir_rad);
    v_out = -speed_ms * std::cos(dir_rad);
}

// Compute mixing ratio from dewpoint temperature using the Tetens formula
static double qv_from_dewpoint(double Td_C, double p_hPa) {
    double es = atm::saturation_vapor_pressure(Td_C); // Pa
    double p_Pa = p_hPa * 100.0;
    double epsilon = atm::Rd / atm::Rv; // ~0.622
    double qv = epsilon * es / (p_Pa - es);
    return std::max(qv, 0.0);
}

// Linear interpolation helper
static double lerp(double t, double a, double b) {
    return a + t * (b - a);
}

void ingest_radiosonde(Grid& grid, const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[radiosonde] Warning: could not open '" << path << "'\n";
        return;
    }

    // Parse CSV
    std::vector<SondeRecord> records;
    std::string line;

    // Skip header line
    if (std::getline(file, line)) {
        // Check if it looks like data (starts with digit or minus)
        if (!line.empty() && (std::isdigit(static_cast<unsigned char>(line[0])) ||
                              line[0] == '-' || line[0] == '+')) {
            std::istringstream iss(line);
            SondeRecord rec;
            char comma;
            if (iss >> rec.p_hPa >> comma >> rec.height_m >> comma
                    >> rec.T_C >> comma >> rec.Td_C >> comma
                    >> rec.wind_dir >> comma >> rec.wind_speed) {
                records.push_back(rec);
            }
        }
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        SondeRecord rec;
        char comma;
        if (iss >> rec.p_hPa >> comma >> rec.height_m >> comma
                >> rec.T_C >> comma >> rec.Td_C >> comma
                >> rec.wind_dir >> comma >> rec.wind_speed) {
            records.push_back(rec);
        } else {
            std::cerr << "[radiosonde] Skipping malformed line: " << line << "\n";
        }
    }

    if (records.empty()) {
        std::cerr << "[radiosonde] No valid records in '" << path << "'\n";
        return;
    }

    std::cout << "[radiosonde] Loaded " << records.size() << " records from '"
              << path << "'\n";

    // Sort by ascending height
    std::sort(records.begin(), records.end(),
              [](const SondeRecord& a, const SondeRecord& b) {
                  return a.height_m < b.height_m;
              });

    // Build interpolated profiles for each grid level
    int center_i = grid.Nx / 2;
    int center_j = grid.Ny / 2;

    // Maximum blending radius: half the smaller horizontal domain extent
    double Lx = grid.Nx * grid.dx;
    double Ly = grid.Ny * grid.dy;
    double blend_radius = 0.4 * std::min(Lx, Ly);

    for (int k = 0; k < grid.Nz; ++k) {
        double z_target = grid.z[k];

        // Interpolate sonde data to this height level
        double T_C_interp = records.front().T_C;
        double Td_C_interp = records.front().Td_C;
        double dir_interp = records.front().wind_dir;
        double spd_interp = records.front().wind_speed;
        double p_hPa_interp = records.front().p_hPa;

        if (z_target <= records.front().height_m) {
            // Below lowest record -- clamp
        } else if (z_target >= records.back().height_m) {
            T_C_interp  = records.back().T_C;
            Td_C_interp = records.back().Td_C;
            dir_interp  = records.back().wind_dir;
            spd_interp  = records.back().wind_speed;
            p_hPa_interp = records.back().p_hPa;
        } else {
            for (std::size_t n = 0; n + 1 < records.size(); ++n) {
                if (z_target >= records[n].height_m &&
                    z_target <= records[n + 1].height_m) {
                    double dz = records[n + 1].height_m - records[n].height_m;
                    double t = (dz > 1e-6) ? (z_target - records[n].height_m) / dz : 0.0;
                    T_C_interp   = lerp(t, records[n].T_C,       records[n + 1].T_C);
                    Td_C_interp  = lerp(t, records[n].Td_C,      records[n + 1].Td_C);
                    dir_interp   = lerp(t, records[n].wind_dir,  records[n + 1].wind_dir);
                    spd_interp   = lerp(t, records[n].wind_speed,records[n + 1].wind_speed);
                    p_hPa_interp = lerp(t, records[n].p_hPa,     records[n + 1].p_hPa);
                    break;
                }
            }
        }

        // Convert sonde values
        double T_K = T_C_interp + 273.15;
        double u_sonde, v_sonde;
        wind_to_uv(dir_interp, spd_interp, u_sonde, v_sonde);
        double qv_sonde = qv_from_dewpoint(Td_C_interp, p_hPa_interp);
        double p_Pa = p_hPa_interp * 100.0;
        double theta_sonde = atm::potential_temperature(T_K, p_Pa);

        // Apply to grid with distance-based blending from the domain center
        double cx = (center_i + 0.5) * grid.dx;
        double cy = (center_j + 0.5) * grid.dy;

        for (int j = 0; j < grid.Ny; ++j) {
            double y = (j + 0.5) * grid.dy;
            for (int i = 0; i < grid.Nx; ++i) {
                double x = (i + 0.5) * grid.dx;
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;

                double dist = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                double weight = 1.0;
                if (dist > blend_radius) {
                    weight = 0.0;
                } else if (dist > 0.5 * blend_radius) {
                    // Smooth cosine taper between 50% and 100% of blend radius
                    double frac = (dist - 0.5 * blend_radius) / (0.5 * blend_radius);
                    weight = 0.5 * (1.0 + std::cos(M_PI * frac));
                }

                if (weight > 0.0) {
                    grid.theta[id] = weight * theta_sonde + (1.0 - weight) * grid.theta[id];
                    grid.qv[id]    = weight * qv_sonde    + (1.0 - weight) * grid.qv[id];
                    grid.u[id]     = weight * u_sonde     + (1.0 - weight) * grid.u[id];
                    grid.v[id]     = weight * v_sonde     + (1.0 - weight) * grid.v[id];
                }
            }
        }
    }

    std::cout << "[radiosonde] Applied radiosonde profile to grid (center blending).\n";
}

} // namespace io
