#include "grid.hpp"

Grid::Grid(const std::array<int,3>& size_) : size(size_) {
    std::size_t N = static_cast<std::size_t>(size[0]) * size[1] * size[2];
    u.assign(N, 0.0);
    v.assign(N, 0.0);
    w.assign(N, 0.0);
}
