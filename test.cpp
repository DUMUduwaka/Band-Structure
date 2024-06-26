#include <iostream>
#include <cmath>
#include <array>
// #include <iomanip>

// Define constants
const double pi = M_PI;
const double theta = 1.2 * pi / 180;         // Twist angle in radians
const double V = 8;                          // Potential coefficient in meV
const double psi = -89.6 * pi / 180;         // Potential angle in radians
const double w = -8.5;                       // Tunneling amplitude in meV
const double Vz = 0;                         // Displacement field in meV
const double meff = 0.62;                    // Electron mass in bare units
const double a0 = 3.472e-10;                 // Lattice constant of MoTe2 in meters
const double area = a0 * a0 * sqrt(3) / 2;   // Area of untwisted case
const double a_M = a0 / theta;               // Moire lattice constant
const double me = 0.51099895e9;              // Mass conversion factor
const double c = 299792458;                  // Speed of light in m/s
const double c_sqd = c * c;                  // Speed of light squared
const double hbar = 4.135667e-12 / (2 * pi); // Planck's reduced constant in meV
const double mstar = meff * me;

const double elec_charge = std::sqrt(1.4399764e9 * 1e-15);
const double e_squared = elec_charge * elec_charge;
const double moire_area = std::sqrt(3) * a_M * a_M / 2;
const int num_k = 50;

const double mu = 22.45;          // Chemical Potential
const double epsilon = 1 / 0.075; // Dielectric constant
const double d_gate = 5e-9;       // Distance to gates for the voltage in meters

const double g_G = 4 * pi / (std::sqrt(3) * a_M); // Factor for G-vectors
const double V_c = 2 * pi * e_squared / (epsilon * g_G);
const double kin = -c_sqd * hbar * hbar * g_G * g_G / (2 * mstar);

int main()
{
    std::cout << "ahdsfs" << std::endl;
    return 0;
}
