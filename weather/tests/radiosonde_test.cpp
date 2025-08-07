#include "../io/radiosonde.hpp"
#include <cassert>
#include <iostream>

int main() {
    auto data = io::ingest_radiosonde("weather/tests/sample_radiosonde.txt");
    assert(data.size() == 2);
    assert(data[0].pressure == 1000);
    assert(data[1].humidity == 60);
    std::cout << "radiosonde test passed" << std::endl;
    return 0;
}
