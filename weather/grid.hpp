#pragma once
#include <array>
#include <vector>

struct Grid {
    std::array<int,3> size;
    std::vector<double> u, v, w;
    explicit Grid(const std::array<int,3>& size);
};
