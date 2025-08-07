#include "config.hpp"
#include "grid.hpp"
#include "terrain.hpp"
#include "physics/thermal.hpp"
#include "physics/cloud.hpp"
#include "io/gfs.hpp"
#include "io/radiosonde.hpp"
#include "io/aws.hpp"
#include "io/output.hpp"
#include <iostream>

int main() {
    Config cfg = load_config("config.toml");
    Grid grid(cfg.domain_size);
    load_terrain(cfg.terrain_file);
    if (cfg.use_gfs) {
        io::load_gfs(cfg.gfs_file);
    }
    io::ingest_radiosonde(cfg.radiosonde_file);
    io::ingest_aws(cfg.aws_file);

    for (int step = 0; step < cfg.run_steps; ++step) {
        physics::apply_thermal_buoyancy();
        physics::update_cloud_microphysics();
        if (step % 100 == 0) {
            io::write_output(cfg.output_dir);
        }
    }
    std::cout << "Simulation complete" << std::endl;
    return 0;
}
