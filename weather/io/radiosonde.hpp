/*
 * Header for radiosonde data ingestion.
 *
 * This exposes a simple `RadiosondeRecord` struct and returns a vector of
 * records when ingesting a radiosonde text file. Each line of the file is
 * expected to contain three whitespace separated values representing the
 * pressure (hPa), temperature (K) and relative humidity (%).
 */

#pragma once

#include <string>
#include <vector>

namespace io {

struct RadiosondeRecord {
    double pressure;
    double temperature;
    double humidity;
};

std::vector<RadiosondeRecord> ingest_radiosonde(const std::string& path);

} // namespace io
