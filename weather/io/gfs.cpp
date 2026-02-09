#include "gfs.hpp"
#include "../constants.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace io {

// A single GFS pressure-level record
struct GFSLevel {
    double p_hPa;      // pressure [hPa]
    double T_K;         // temperature [K]
    double u_ms;        // u-wind [m/s]
    double v_ms;        // v-wind [m/s]
    double qv_kgkg;     // specific humidity [kg/kg]
};

// Approximate geopotential height from pressure using the hypsometric equation
// relative to the standard atmosphere.
static double pressure_to_height(double p_hPa) {
    // H = (Rd * T0 / g) * ln(p0 / p)  -- isothermal approximation
    double p_Pa = p_hPa * 100.0;
    if (p_Pa <= 0.0) p_Pa = 1.0;
    double H = (atm::Rd * atm::T0 / atm::g) * std::log(atm::p0 / p_Pa);
    return std::max(H, 0.0);
}

// Log-pressure linear interpolation between two levels
static double interp_logp(double z_target,
                          double z0, double z1,
                          double val0, double val1) {
    if (std::abs(z1 - z0) < 1e-6) return val0;
    double t = (z_target - z0) / (z1 - z0);
    t = std::max(0.0, std::min(1.0, t));
    return val0 + t * (val1 - val0);
}

// Initialise the grid with standard-atmosphere defaults when no file is available
static void apply_standard_atmosphere(Grid& grid) {
    std::cout << "[gfs] Applying standard atmosphere defaults.\n";
    for (int k = 0; k < grid.Nz; ++k) {
        double z_m = grid.z[k];
        double T_k = atm::T0 - atm::gamma_d * z_m;
        if (T_k < 200.0) T_k = 200.0; // clamp at stratosphere floor
        double p_k = atm::p0 * std::pow(T_k / atm::T0, atm::g / (atm::Rd * atm::gamma_d));
        double theta_k = atm::potential_temperature(T_k, p_k);
        double qv_k = 0.005 * std::exp(-z_m / 2500.0); // simple moisture decrease

        for (int j = 0; j < grid.Ny; ++j) {
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;
                grid.theta[id] = theta_k;
                grid.qv[id]    = std::max(qv_k, 0.0);
                grid.u[id]     = 0.0;
                grid.v[id]     = 0.0;
                grid.w[id]     = 0.0;
                grid.p[id]     = p_k;
                grid.T[id]     = T_k;
                grid.rho[id]   = atm::density(p_k, T_k);
            }
        }
    }
}

void load_gfs(Grid& grid, const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[gfs] Warning: could not open '" << path
                  << "'. Using standard atmosphere defaults.\n";
        apply_standard_atmosphere(grid);
        return;
    }

    // Parse CSV records
    std::vector<GFSLevel> levels;
    std::string line;

    // Skip header line
    if (std::getline(file, line)) {
        // check if first line looks like data (starts with a digit)
        if (!line.empty() && (std::isdigit(static_cast<unsigned char>(line[0])) ||
                              line[0] == '-' || line[0] == '+')) {
            // No header -- parse this line as data
            std::istringstream iss(line);
            GFSLevel lev;
            char comma;
            if (iss >> lev.p_hPa >> comma >> lev.T_K >> comma
                    >> lev.u_ms  >> comma >> lev.v_ms >> comma >> lev.qv_kgkg) {
                levels.push_back(lev);
            }
        }
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        GFSLevel lev;
        char comma;
        if (iss >> lev.p_hPa >> comma >> lev.T_K >> comma
                >> lev.u_ms  >> comma >> lev.v_ms >> comma >> lev.qv_kgkg) {
            levels.push_back(lev);
        } else {
            std::cerr << "[gfs] Skipping malformed line: " << line << "\n";
        }
    }

    if (levels.empty()) {
        std::cerr << "[gfs] No valid records found in '" << path
                  << "'. Using standard atmosphere defaults.\n";
        apply_standard_atmosphere(grid);
        return;
    }

    std::cout << "[gfs] Loaded " << levels.size() << " pressure levels from '"
              << path << "'\n";

    // Compute approximate height for each GFS level
    struct HeightRecord {
        double z_m, T_K, u, v, qv;
    };
    std::vector<HeightRecord> profile(levels.size());
    for (std::size_t n = 0; n < levels.size(); ++n) {
        profile[n].z_m = pressure_to_height(levels[n].p_hPa);
        profile[n].T_K = levels[n].T_K;
        profile[n].u   = levels[n].u_ms;
        profile[n].v   = levels[n].v_ms;
        profile[n].qv  = levels[n].qv_kgkg;
    }

    // Sort profile by ascending height
    std::sort(profile.begin(), profile.end(),
              [](const HeightRecord& a, const HeightRecord& b) {
                  return a.z_m < b.z_m;
              });

    // Interpolate to each grid level and apply uniformly across the domain
    for (int k = 0; k < grid.Nz; ++k) {
        double z_target = grid.z[k];

        double T_interp  = profile.front().T_K;
        double u_interp  = profile.front().u;
        double v_interp  = profile.front().v;
        double qv_interp = profile.front().qv;

        if (z_target <= profile.front().z_m) {
            // Below the lowest level -- clamp
            T_interp  = profile.front().T_K;
            u_interp  = profile.front().u;
            v_interp  = profile.front().v;
            qv_interp = profile.front().qv;
        } else if (z_target >= profile.back().z_m) {
            // Above the highest level -- clamp
            T_interp  = profile.back().T_K;
            u_interp  = profile.back().u;
            v_interp  = profile.back().v;
            qv_interp = profile.back().qv;
        } else {
            // Find bracketing levels
            for (std::size_t n = 0; n + 1 < profile.size(); ++n) {
                if (z_target >= profile[n].z_m && z_target <= profile[n + 1].z_m) {
                    T_interp  = interp_logp(z_target, profile[n].z_m,
                                            profile[n + 1].z_m,
                                            profile[n].T_K, profile[n + 1].T_K);
                    u_interp  = interp_logp(z_target, profile[n].z_m,
                                            profile[n + 1].z_m,
                                            profile[n].u, profile[n + 1].u);
                    v_interp  = interp_logp(z_target, profile[n].z_m,
                                            profile[n + 1].z_m,
                                            profile[n].v, profile[n + 1].v);
                    qv_interp = interp_logp(z_target, profile[n].z_m,
                                            profile[n + 1].z_m,
                                            profile[n].qv, profile[n + 1].qv);
                    break;
                }
            }
        }

        // Estimate pressure at this height (standard atmosphere relation)
        double p_est = atm::p0 * std::pow(
            std::max(1e-6, 1.0 - atm::gamma_d * z_target / atm::T0),
            atm::g / (atm::Rd * atm::gamma_d));
        double theta_val = atm::potential_temperature(T_interp, p_est);

        // Apply uniformly across horizontal domain
        for (int j = 0; j < grid.Ny; ++j) {
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;
                grid.theta[id] = theta_val;
                grid.qv[id]    = std::max(qv_interp, 0.0);
                grid.u[id]     = u_interp;
                grid.v[id]     = v_interp;
                grid.w[id]     = 0.0;
                grid.p[id]     = p_est;
                grid.T[id]     = T_interp;
                grid.rho[id]   = atm::density(p_est, T_interp);
            }
        }
    }

    std::cout << "[gfs] Applied GFS profile to grid.\n";
}

} // namespace io
