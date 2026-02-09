#include "output.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cstring>

namespace io {

// ---------------------------------------------------------------------------
// ensure_directory -- create directory tree (equivalent of mkdir -p)
// ---------------------------------------------------------------------------
void ensure_directory(const std::string& dir) {
    if (dir.empty()) return;

    // Walk through the path components and create each one
    std::string accumulated;
    for (std::size_t i = 0; i < dir.size(); ++i) {
        accumulated += dir[i];
        if (dir[i] == '/' || i == dir.size() - 1) {
            struct stat st;
            if (stat(accumulated.c_str(), &st) != 0) {
                if (mkdir(accumulated.c_str(), 0755) != 0 && errno != EEXIST) {
                    std::cerr << "[output] Failed to create directory '"
                              << accumulated << "': " << std::strerror(errno) << "\n";
                    return;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// write_vtk -- VTK legacy structured-points file
// ---------------------------------------------------------------------------
void write_vtk(const Grid& grid, const std::string& dir, int step) {
    ensure_directory(dir);

    std::ostringstream fname;
    fname << dir << "/weather_" << std::setfill('0') << std::setw(6) << step << ".vtk";

    std::ofstream out(fname.str());
    if (!out.is_open()) {
        std::cerr << "[output] Could not open VTK file '" << fname.str() << "'\n";
        return;
    }

    // Use the average dz for spacing (grid has variable dz, but VTK
    // STRUCTURED_POINTS requires uniform spacing)
    double avg_dz = 0.0;
    for (int k = 0; k < grid.Nz; ++k) avg_dz += grid.dz[k];
    avg_dz /= grid.Nz;

    // Header
    out << "# vtk DataFile Version 3.0\n";
    out << "Weather Step " << step << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << grid.Nx << " " << grid.Ny << " " << grid.Nz << "\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << grid.dx << " " << grid.dy << " " << avg_dz << "\n";
    out << "POINT_DATA " << grid.total() << "\n";

    // Temperature
    out << "SCALARS temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int k = 0; k < grid.Nz; ++k)
        for (int j = 0; j < grid.Ny; ++j)
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                out << (grid.solid[id] ? 0.0 : grid.T[id]) << "\n";
            }

    // Wind speed (scalar magnitude)
    out << "SCALARS wind_speed double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int k = 0; k < grid.Nz; ++k)
        for (int j = 0; j < grid.Ny; ++j)
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) {
                    out << 0.0 << "\n";
                } else {
                    double spd = std::sqrt(grid.u[id] * grid.u[id] +
                                           grid.v[id] * grid.v[id] +
                                           grid.w[id] * grid.w[id]);
                    out << spd << "\n";
                }
            }

    // Pressure
    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int k = 0; k < grid.Nz; ++k)
        for (int j = 0; j < grid.Ny; ++j)
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                out << (grid.solid[id] ? 0.0 : grid.p[id]) << "\n";
            }

    // Cloud water
    out << "SCALARS cloud_water double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int k = 0; k < grid.Nz; ++k)
        for (int j = 0; j < grid.Ny; ++j)
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                out << (grid.solid[id] ? 0.0 : grid.qc[id]) << "\n";
            }

    // Rain water
    out << "SCALARS rain_water double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int k = 0; k < grid.Nz; ++k)
        for (int j = 0; j < grid.Ny; ++j)
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                out << (grid.solid[id] ? 0.0 : grid.qr[id]) << "\n";
            }

    // Velocity (vector field)
    out << "VECTORS velocity double\n";
    for (int k = 0; k < grid.Nz; ++k)
        for (int j = 0; j < grid.Ny; ++j)
            for (int i = 0; i < grid.Nx; ++i) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) {
                    out << "0 0 0\n";
                } else {
                    out << grid.u[id] << " " << grid.v[id] << " " << grid.w[id] << "\n";
                }
            }

    out.close();
    std::cout << "[output] Wrote VTK file: " << fname.str() << "\n";
}

// ---------------------------------------------------------------------------
// write_csv_summary -- subsampled CSV snapshot
// ---------------------------------------------------------------------------
void write_csv_summary(const Grid& grid, const std::string& dir, int step) {
    ensure_directory(dir);

    std::ostringstream fname;
    fname << dir << "/summary_" << std::setfill('0') << std::setw(6) << step << ".csv";

    std::ofstream out(fname.str());
    if (!out.is_open()) {
        std::cerr << "[output] Could not open CSV file '" << fname.str() << "'\n";
        return;
    }

    // Header
    out << "x,y,z,T,p,u,v,w,qv,qc,qr,theta\n";
    out << std::scientific << std::setprecision(6);

    // Subsample: every 4th cell in each direction
    const int stride = 4;
    for (int k = 0; k < grid.Nz; k += stride) {
        for (int j = 0; j < grid.Ny; j += stride) {
            for (int i = 0; i < grid.Nx; i += stride) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.solid[id]) continue;

                double x = (i + 0.5) * grid.dx;
                double y = (j + 0.5) * grid.dy;
                double z = grid.z[k];

                out << x << "," << y << "," << z << ","
                    << grid.T[id]     << "," << grid.p[id]  << ","
                    << grid.u[id]     << "," << grid.v[id]  << "," << grid.w[id] << ","
                    << grid.qv[id]    << "," << grid.qc[id] << "," << grid.qr[id] << ","
                    << grid.theta[id] << "\n";
            }
        }
    }

    out.close();
    std::cout << "[output] Wrote CSV summary: " << fname.str() << "\n";
}

} // namespace io
