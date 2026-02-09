#include "terrain.hpp"
#include "constants.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <vector>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Read a PGM image into a flat vector of pixel values (0..maxval)
static bool read_pgm(const std::string& path, int& img_w, int& img_h,
                     int& maxval, std::vector<int>& pixels) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return false;

    std::string magic;
    file >> magic;
    if (magic != "P2" && magic != "P5") {
        std::cerr << "[terrain] Unsupported PGM magic: " << magic << "\n";
        return false;
    }

    // Skip comments
    auto skip_comments = [&]() {
        while (file.peek() == '#' || file.peek() == '\n' || file.peek() == '\r' ||
               file.peek() == ' ') {
            if (file.peek() == '#') {
                std::string line;
                std::getline(file, line);
            } else {
                file.get();
            }
        }
    };

    skip_comments();
    file >> img_w;
    skip_comments();
    file >> img_h;
    skip_comments();
    file >> maxval;

    if (img_w <= 0 || img_h <= 0 || maxval <= 0) {
        std::cerr << "[terrain] Invalid PGM dimensions: " << img_w << "x" << img_h
                  << " maxval=" << maxval << "\n";
        return false;
    }

    pixels.resize((std::size_t)img_w * img_h);

    if (magic == "P2") {
        // ASCII format
        for (std::size_t i = 0; i < pixels.size(); ++i) {
            if (!(file >> pixels[i])) {
                std::cerr << "[terrain] Premature end of P2 data at pixel " << i << "\n";
                return false;
            }
        }
    } else {
        // P5 binary format
        file.get(); // consume single whitespace after maxval
        if (maxval <= 255) {
            for (std::size_t i = 0; i < pixels.size(); ++i) {
                uint8_t byte;
                file.read(reinterpret_cast<char*>(&byte), 1);
                pixels[i] = byte;
            }
        } else {
            // 16-bit big-endian
            for (std::size_t i = 0; i < pixels.size(); ++i) {
                uint8_t hi, lo;
                file.read(reinterpret_cast<char*>(&hi), 1);
                file.read(reinterpret_cast<char*>(&lo), 1);
                pixels[i] = (hi << 8) | lo;
            }
        }
        if (!file) {
            std::cerr << "[terrain] Premature end of P5 data\n";
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

void load_terrain(Grid& grid, const std::string& path) {
    int img_w = 0, img_h = 0, maxval = 0;
    std::vector<int> pixels;

    if (!read_pgm(path, img_w, img_h, maxval, pixels)) {
        std::cerr << "[terrain] Could not load PGM file '" << path
                  << "', falling back to procedural terrain.\n";
        generate_terrain(grid, 1500.0);
        return;
    }

    std::cout << "[terrain] Loaded PGM " << img_w << "x" << img_h
              << " maxval=" << maxval << " from '" << path << "'\n";

    const double max_terrain_height = 3000.0; // metres
    const double scale = max_terrain_height / static_cast<double>(maxval);

    // Resample to Nx * Ny using nearest-neighbour interpolation
    for (int j = 0; j < grid.Ny; ++j) {
        for (int i = 0; i < grid.Nx; ++i) {
            int src_x = static_cast<int>((static_cast<double>(i) / grid.Nx) * img_w);
            int src_y = static_cast<int>((static_cast<double>(j) / grid.Ny) * img_h);
            src_x = std::min(src_x, img_w - 1);
            src_y = std::min(src_y, img_h - 1);

            double h = pixels[(std::size_t)src_y * img_w + src_x] * scale;
            grid.terrain_height[(std::size_t)j * grid.Nx + i] = h;
        }
    }

    apply_terrain_to_grid(grid);
}

void generate_terrain(Grid& grid, double max_height) {
    const int num_hills = 4;
    // Reproducible pseudo-random placement via fixed seed arithmetic
    const unsigned seed = 42u;

    struct Hill {
        double cx, cy, sigma, height;
    };

    // Domain physical extents
    double Lx = grid.Nx * grid.dx;
    double Ly = grid.Ny * grid.dy;

    // Generate hills deterministically from the seed
    std::vector<Hill> hills(num_hills);
    unsigned rng = seed;
    auto next_rand = [&]() -> double {
        // Simple LCG producing values in [0,1)
        rng = rng * 1103515245u + 12345u;
        return static_cast<double>((rng >> 16) & 0x7FFF) / 32768.0;
    };

    for (int h = 0; h < num_hills; ++h) {
        hills[h].cx     = 0.1 * Lx + 0.8 * Lx * next_rand();
        hills[h].cy     = 0.1 * Ly + 0.8 * Ly * next_rand();
        hills[h].sigma  = 0.05 * std::min(Lx, Ly) + 0.15 * std::min(Lx, Ly) * next_rand();
        hills[h].height = 0.3 * max_height + 0.7 * max_height * next_rand();
    }

    // Evaluate the terrain height field
    for (int j = 0; j < grid.Ny; ++j) {
        double y = (j + 0.5) * grid.dy;
        for (int i = 0; i < grid.Nx; ++i) {
            double x = (i + 0.5) * grid.dx;
            double elevation = 0.0;
            for (int h = 0; h < num_hills; ++h) {
                double ddx = x - hills[h].cx;
                double ddy = y - hills[h].cy;
                double r2  = ddx * ddx + ddy * ddy;
                double s2  = 2.0 * hills[h].sigma * hills[h].sigma;
                elevation += hills[h].height * std::exp(-r2 / s2);
            }
            grid.terrain_height[(std::size_t)j * grid.Nx + i] = elevation;
        }
    }

    std::cout << "[terrain] Generated procedural terrain with " << num_hills
              << " Gaussian hills, max_height=" << max_height << " m\n";

    apply_terrain_to_grid(grid);
}

void apply_terrain_to_grid(Grid& grid) {
    int solid_count = 0;
    for (int j = 0; j < grid.Ny; ++j) {
        for (int i = 0; i < grid.Nx; ++i) {
            double th = grid.terrain_height[(std::size_t)j * grid.Nx + i];
            for (int k = 0; k < grid.Nz; ++k) {
                std::size_t id = grid.idx(i, j, k);
                if (grid.z[k] < th) {
                    grid.solid[id] = true;
                    grid.u[id] = 0.0;
                    grid.v[id] = 0.0;
                    grid.w[id] = 0.0;
                    ++solid_count;
                } else {
                    grid.solid[id] = false;
                }
            }
        }
    }

    std::cout << "[terrain] Applied terrain to grid: " << solid_count << " / "
              << grid.total() << " cells marked solid ("
              << (100.0 * solid_count / grid.total()) << "%)\n";
}
