// main.cpp
#include <iostream>
#include <vector>
#include <array>
#include <complex>
#include "constants.h"
#include "utilities.h"

// first shell is just the 0 vector
double g_0[2] = {0 * K_1[0] + 0 * K_2[0], 0 * K_1[1] + 0 * K_2[1]};

// second shell is the first six vectors of equal magnitude
double g_1[2] = {0 * K_1[0] + 1 * K_2[0], 0 * K_1[1] + 1 * K_2[1]};
double g_2[2] = {1 * K_1[0] + 0 * K_2[0], 1 * K_1[1] + 0 * K_2[1]};
double g_3[2] = {1 * K_1[0] + -1 * K_2[0], 1 * K_1[1] + -1 * K_2[1]};
double g_4[2] = {0 * K_1[0] + -1 * K_2[0], 0 * K_1[1] + -1 * K_2[1]};
double g_5[2] = {-1 * K_1[0] + 0 * K_2[0], -1 * K_1[1] + 0 * K_2[1]};
double g_6[2] = {-1 * K_1[0] + 1 * K_2[0], -1 * K_1[1] + 1 * K_2[1]};

int main()
{
    std::cout << "jhdh" << std::endl;
    for (int i = 0; i < 2; i++)
    {
        std::cout << g_6[i] << " ";
    }
    std::cout << std::endl;
    return 0;
}