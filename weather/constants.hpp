#pragma once
#include <cmath>

namespace atm {

// Fundamental physical constants
constexpr double g       = 9.80665;       // gravitational acceleration [m/s^2]
constexpr double Rd      = 287.05;        // gas constant for dry air [J/(kg*K)]
constexpr double Rv      = 461.5;         // gas constant for water vapor [J/(kg*K)]
constexpr double cp_d    = 1004.0;        // specific heat of dry air at const pressure [J/(kg*K)]
constexpr double cv_d    = 717.0;         // specific heat of dry air at const volume [J/(kg*K)]
constexpr double cp_v    = 1870.0;        // specific heat of water vapor at const pressure [J/(kg*K)]
constexpr double Lv      = 2.501e6;       // latent heat of vaporization at 0C [J/kg]
constexpr double Ls      = 2.834e6;       // latent heat of sublimation [J/kg]
constexpr double Lf      = 3.34e5;        // latent heat of fusion [J/kg]
constexpr double Omega   = 7.2921e-5;     // Earth's angular velocity [rad/s]
constexpr double R_earth = 6.371e6;       // Earth's mean radius [m]
constexpr double p0      = 100000.0;      // reference pressure [Pa]
constexpr double T0      = 288.15;        // standard sea-level temperature [K]
constexpr double rho0    = 1.225;         // standard sea-level density [kg/m^3]
constexpr double kappa   = Rd / cp_d;     // Poisson constant (R/cp) ~ 0.2857

// Dry adiabatic lapse rate [K/m]
constexpr double gamma_d = g / cp_d;      // ~0.00976 K/m

// Saturation vapor pressure constants (Tetens formula)
constexpr double es0     = 611.2;         // reference saturation vapor pressure at 0C [Pa]
constexpr double Tes_a   = 17.67;         // Tetens coefficient a
constexpr double Tes_b   = 243.5;         // Tetens coefficient b [C]

// Turbulence constants
constexpr double Km_default = 500.0;      // default horizontal momentum diffusion [m^2/s]
constexpr double Kh_default = 500.0;      // default horizontal thermal diffusion [m^2/s]
constexpr double Cs         = 0.2;        // Smagorinsky constant

// Kessler microphysics constants
constexpr double qc_autoconv_threshold = 1.0e-3; // autoconversion threshold [kg/kg]
constexpr double autoconv_rate         = 1.0e-3; // autoconversion rate [1/s]
constexpr double accretion_rate        = 2.2;     // accretion rate constant
constexpr double evaporation_rate      = 1.0e-3;  // rain evaporation rate
constexpr double terminal_velocity_coeff = 5.0;   // Vt = coeff * (rho0/rho)^0.5 * qr^0.1375

// Radiation constants
constexpr double stefan_boltzmann = 5.67e-8; // Stefan-Boltzmann constant [W/(m^2*K^4)]
constexpr double solar_constant  = 1361.0;   // solar constant [W/m^2]

// Saturation vapor pressure over water (Tetens formula)
inline double saturation_vapor_pressure(double T_celsius) {
    return es0 * std::exp(Tes_a * T_celsius / (T_celsius + Tes_b));
}

// Saturation mixing ratio
inline double saturation_mixing_ratio(double p_Pa, double T_K) {
    double T_C = T_K - 273.15;
    double es = saturation_vapor_pressure(T_C);
    double epsilon = Rd / Rv; // ~0.622
    return epsilon * es / (p_Pa - es);
}

// Potential temperature from temperature and pressure
inline double potential_temperature(double T_K, double p_Pa) {
    return T_K * std::pow(p0 / p_Pa, kappa);
}

// Temperature from potential temperature and pressure
inline double temperature_from_theta(double theta, double p_Pa) {
    return theta * std::pow(p_Pa / p0, kappa);
}

// Virtual temperature (accounting for moisture)
inline double virtual_temperature(double T_K, double qv) {
    return T_K * (1.0 + 0.608 * qv);
}

// Coriolis parameter at given latitude [rad]
inline double coriolis_parameter(double latitude_rad) {
    return 2.0 * Omega * std::sin(latitude_rad);
}

// Dry air density from ideal gas law
inline double density(double p_Pa, double T_K) {
    return p_Pa / (Rd * T_K);
}

// Pressure from density and temperature
inline double pressure(double rho_val, double T_K) {
    return rho_val * Rd * T_K;
}

// Moist static energy
inline double moist_static_energy(double T_K, double z, double qv) {
    return cp_d * T_K + g * z + Lv * qv;
}

} // namespace atm
