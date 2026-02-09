#include "../config.hpp"
#include <cassert>
#include <iostream>
#include <cmath>

int main() {
    Config cfg = load_config("weather/config.toml");

    assert(cfg.nx == 80);
    assert(cfg.ny == 80);
    assert(cfg.nz == 40);
    assert(std::abs(cfg.dx - 2000.0) < 0.01);
    assert(std::abs(cfg.dy - 2000.0) < 0.01);
    assert(std::abs(cfg.dz_base - 250.0) < 0.01);
    assert(std::abs(cfg.timestep - 2.0) < 0.01);
    assert(cfg.run_steps == 1800);
    assert(cfg.output_interval == 60);
    assert(std::abs(cfg.latitude - (-27.5)) < 0.01);
    assert(!cfg.use_gfs);
    assert(cfg.add_warm_bubble);
    assert(std::abs(cfg.bubble_dtheta - 3.0) < 0.01);
    assert(cfg.output_vtk);
    assert(cfg.output_csv);

    std::cout << "config test passed" << std::endl;
    return 0;
}
