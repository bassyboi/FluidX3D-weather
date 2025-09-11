#include "grid.hpp"
#include <limits>
#include <stdexcept>

Grid::Grid(const std::array<int,3>& size_) : size(size_) {
    // Check for overflow before multiplying
    std::size_t s0 = static_cast<std::size_t>(size[0]);
    std::size_t s1 = static_cast<std::size_t>(size[1]);
    std::size_t s2 = static_cast<std::size_t>(size[2]);
    std::size_t max = std::numeric_limits<std::size_t>::max();
    if (s1 != 0 && s0 > max / s1) {
        throw std::overflow_error("Grid dimension overflow (size[0] * size[1])");
    }
    std::size_t prod01 = s0 * s1;
    if (s2 != 0 && prod01 > max / s2) {
        throw std::overflow_error("Grid dimension overflow (size[0] * size[1] * size[2])");
    }
    std::size_t N = prod01 * s2;
    u.assign(N, 0.0);
    v.assign(N, 0.0);
    w.assign(N, 0.0);
}
