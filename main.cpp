// main.cpp
#include <iostream>
#include <vector>
#include <array>
#include <complex>
#include "constants.h"
#include <armadillo>

using namespace std;
using namespace arma;

// Reciprocal lattice vectors
rowvec K_1 = g_G * rowvec({0.5, std::sqrt(3) / 2});
rowvec K_2 = g_G * rowvec({1.0, 0.0});

mat ns = {
    {1, 0},
    {0, 1},
    {-1, 1},
    {-1, 0},
    {0, -1},
    {1, -1}};
// first shell is just the 0 vector
rowvec g_0 = 0 * K_1 + 0 * K_2;

// second shell is the first six vectors of equal magnitude
rowvec g_1 = 0 * K_1 + 1 * K_2;
rowvec g_2 = 1 * K_1 + 0 * K_2;
rowvec g_3 = 1 * K_1 + -1 * K_2;
rowvec g_4 = 0 * K_1 + -1 * K_2;
rowvec g_5 = -1 * K_1 + 0 * K_2;
rowvec g_6 = -1 * K_1 + 1 * K_2;

// third shell can be constructed from specific sums of the first
rowvec g_7 = g_1 + g_2;
rowvec g_8 = g_2 + g_3;
rowvec g_9 = g_3 + g_4;
rowvec g_10 = g_4 + g_5;
rowvec g_11 = g_5 + g_6;
rowvec g_12 = g_6 + g_1;

// fourth:
rowvec g_13 = 2 * g_1;
rowvec g_14 = 2 * g_2;
rowvec g_15 = 2 * g_3;
rowvec g_16 = 2 * g_4;
rowvec g_17 = 2 * g_5;
rowvec g_18 = 2 * g_6;

// fifth:
rowvec g_19 = g_7 + g_1;
rowvec g_20 = g_7 + g_2;

rowvec g_21 = g_8 + g_2;
rowvec g_22 = g_8 + g_3;

rowvec g_23 = g_9 + g_3;
rowvec g_24 = g_9 + g_4;

rowvec g_25 = g_10 + g_4;
rowvec g_26 = g_10 + g_5;

rowvec g_27 = g_11 + g_5;
rowvec g_28 = g_11 + g_6;

rowvec g_29 = g_12 + g_6;
rowvec g_30 = g_12 + g_1;

// sixth:
rowvec g_31 = 3 * g_1;
rowvec g_32 = 3 * g_2;
rowvec g_33 = 3 * g_3;
rowvec g_34 = 3 * g_4;
rowvec g_35 = 3 * g_5;
rowvec g_36 = 3 * g_6;

// seventh
rowvec g_37 = 2 * g_7;
rowvec g_38 = 2 * g_8;
rowvec g_39 = 2 * g_9;
rowvec g_40 = 2 * g_10;
rowvec g_41 = 2 * g_11;
rowvec g_42 = 2 * g_12;

std::vector<arma::rowvec> second_shell_gs = {g_0, g_1, g_2, g_3, g_4, g_5, g_6};
std::vector<arma::rowvec> third_shell_gs = {g_0, g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9, g_10, g_11, g_12};
std::vector<arma::rowvec> fourth_shell_gs = {g_0, g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9, g_10, g_11, g_12, g_13, g_14, g_15, g_16, g_17, g_18};
std::vector<arma::rowvec> sixth_shell_gs = {g_0, g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9, g_10, g_11, g_12, g_13, g_14, g_15, g_16, g_17, g_18,
                                            g_19, g_20, g_21, g_22, g_23, g_24, g_25, g_26, g_27, g_28, g_29, g_30, g_31, g_32, g_33, g_34, g_35, g_36};
std::vector<arma::rowvec> seventh_shell_gs = {g_0, g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9, g_10, g_11, g_12, g_13, g_14, g_15, g_16, g_17, g_18,
                                              g_19, g_20, g_21, g_22, g_23, g_24, g_25, g_26, g_27, g_28, g_29, g_30, g_31, g_32, g_33, g_34, g_35,
                                              g_36, g_37, g_38, g_39, g_40, g_41, g_42};

int main()
{
    std::cout << "jhdh" << std::endl;
    for (int i = 0; i < 12; i++)
    {
        std::cout << third_shell_gs[i] << " ";
    }
    std::cout << std::endl;
    return 0;
}
