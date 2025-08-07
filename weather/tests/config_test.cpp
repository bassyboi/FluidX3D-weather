#include "../config.hpp"
#include <cassert>
#include <iostream>

int main() {
    Config cfg = load_config("weather/config.toml");
    assert(cfg.domain_size[0] == 100);
    assert(cfg.domain_size[2] == 60);
    assert(cfg.timestep == 1.0);
    assert(cfg.use_gfs);
    assert(cfg.run_steps == 3600);
    std::cout << "config test passed" << std::endl;
    return 0;
}
