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

// Calculate g_index
arma::vec g_index = arma::regspace(0, sixth_shell_gs.size() - 1);

const int nvalleys = 2;
const int nlayers = 2;
const int ngs = sixth_shell_gs.size();
const int Hdim = nvalleys * nvalleys * ngs;
const int tau = 1;

arma::vec valleys = arma::regspace(0, nvalleys - 1);
arma::vec gs = arma::regspace(0, ngs - 1);
arma::vec layers = arma::regspace(0, nlayers - 1);

// Defining high-symmetric points of the moire Brillouin zone
rowvec Kappa = (4 * pi / (3 * a_M)) * rowvec({1, 0});
rowvec Gamma = rowvec({0, 0});
rowvec Km = rowvec({g_4(0) / 2, -g_4(0) * tan(30 * pi / 180) / 2});
rowvec Kp = rowvec({g_4(0) / 2, g_4(0) * tan(30 * pi / 180) / 2});

// Kinetic Energy Function Definition
cx_mat KE_SH(double kx, double ky, const mat &G)
{
    vec ke_along_g;
    ke_along_g.set_size(2 * G.n_rows); // Ensure size is double the number of rows in G

    cx_mat LAYER_KE_G = eye<cx_mat>(2 * G.n_rows, 2 * G.n_rows);

    for (size_t g = 0; g < G.n_rows; ++g)
    {
        rowvec ke_g_b_vec = {kx + G(g, 0) / g_G - Kp(0) / g_G, ky + G(g, 1) / g_G - Kp(1) / g_G};
        double ke_g_b = dot(ke_g_b_vec, ke_g_b_vec);

        rowvec ke_g_t_vec = {kx + G(g, 0) / g_G - Km(0) / g_G, ky + G(g, 1) / g_G - Km(1) / g_G};
        double ke_g_t = dot(ke_g_t_vec, ke_g_t_vec);

        ke_along_g(2 * g) = ke_g_b;
        ke_along_g(2 * g + 1) = ke_g_t;
    }

    arma::cx_mat ke_along_mat(2 * G.n_rows, 2 * G.n_rows, arma::fill::zeros);
    for (int i = 0; i < 2 * G.n_rows; i++)
    {
        ke_along_mat(i, i) = ke_along_g(i);
    }
    // Convert ke_along_g to a column vector for matrix multiplication
    // cx_vec ke_along_g_vec(ke_along_g);

    return kin * ke_along_mat;
}
int main()
{
    // First Path
    double S1_x = -Km(0);
    double S1_y = -Km(1);
    double F1_x = Gamma(0);
    double F1_y = Gamma(1);

    arma::vec P1_x = arma::linspace(S1_x, F1_x, num_k);
    arma::vec P1_y = arma::zeros(P1_x.size());

    double x = 0;
    double y = 0;
    for (size_t j = 0; j < P1_x.n_elem; j++)
    {
        x = P1_x(j);
        y = x * (F1_y - S1_y) / (F1_x - S1_x);
        P1_y(j) = y;
    }
    // End First Path

    // Second Path
    double S2_x = Gamma(0);
    double S2_y = Gamma(1);
    double F2_x = Km(0);
    double F2_y = Km(1);

    arma::vec P2_x = arma::linspace(S2_x, F2_x, num_k);
    arma::vec P2_y = arma::zeros(P2_x.size());

    double x1 = 0;
    double y1 = 0;
    for (size_t j = 0; j < P2_x.n_elem; j++)
    {
        x1 = P2_x(j);
        y1 = x1 * (F2_y - S2_y) / (F2_x - S2_x);
        P2_y(j) = y1;
    }
    // End of the Second Path

    // Third Path
    double S3_x = Km(0);
    double S3_y = Km(1);
    double F3_x = Kp(0);
    double F3_y = Kp(1);

    arma::vec P3_x = arma::linspace(S3_x, F3_x, num_k);
    arma::vec P3_y = arma::linspace(S3_y, F3_y, num_k);
    // End of the third path

    // Fourth Path
    double S4_x = Kp(0);
    double S4_y = Kp(1);
    double F4_x = -Km(0);
    double F4_y = -Km(1);

    arma::vec P4_x = arma::linspace(S4_x, F4_x, num_k);
    arma::vec P4_y = S4_y * arma::ones(P4_x.size());
    // End of the Fourth path

    for (int i = 0; i < P4_y.size(); i++)
    {
        std::cout << P4_y(i) << std::endl;
    }
    // Print the results
    // P4_x.print("P2_y:");

    const double tol = 1e-5;
    // arma::mat w_on_g2 = arma::zeros()
    size_t rows = sixth_shell_gs.size();
    size_t cols = sixth_shell_gs.size();

    // Create complex matrices
    arma::cx_mat w_on_g2(rows, cols, arma::fill::zeros);
    arma::cx_mat w_on_g3(rows, cols, arma::fill::zeros);
    arma::cx_mat Delta_Odd(rows, cols, arma::fill::zeros);
    arma::cx_mat Delta_Even(rows, cols, arma::fill::zeros);

    w_on_g2.print("w_on_g2:");
    // Print dimensions to verify
    std::cout << "Dimensions: " << rows << "x" << cols << std::endl;

    arma::mat sixth_shell_mat(sixth_shell_gs.size(), g_0.size());
    for (int i = 0; i < sixth_shell_gs.size(); i++)
    {
        sixth_shell_mat.row(i) = sixth_shell_gs[i];
    }
    sixth_shell_mat.print();

    // rowvec g = sixth_shell_mat.row(1);
    // rowvec gp = sixth_shell_mat.row(1);
    // rowvec g_test = g - gp;
    // g_test.print();

    for (size_t G = 0; G < rows; ++G)
    {
        for (size_t Gp = 0; Gp < rows; ++Gp)
        {
            rowvec g = sixth_shell_mat.row(G);
            rowvec gp = sixth_shell_mat.row(Gp);
            rowvec g_test = g - gp;

            if (approx_equal(g_test, g_1, "absdiff", tol) || approx_equal(g_test, g_3, "absdiff", tol) || approx_equal(g_test, g_5, "absdiff", tol))
            {
                Delta_Odd(G, Gp) = cx_double(1, 0);
            }
            if (approx_equal(g_test, g_2, "absdiff", tol) || approx_equal(g_test, g_4, "absdiff", tol) || approx_equal(g_test, g_6, "absdiff", tol))
            {
                Delta_Even(G, Gp) = cx_double(1, 0);
            }
            if (approx_equal(g_test, g_2, "absdiff", tol))
            {
                w_on_g2(G, Gp) = cx_double(1, 0);
            }
            if (approx_equal(g_test, g_3, "absdiff", tol))
            {
                w_on_g3(G, Gp) = cx_double(1, 0);
            }
        }
    }

    // LP matrix
    arma::cx_mat LP(2, 2);
    LP(0, 0) = V * exp(complex<double>(0.0, psi));
    LP(0, 1) = complex<double>(0.0, 0.0);
    LP(1, 0) = complex<double>(0.0, 0.0);
    LP(1, 1) = V * exp(complex<double>(0.0, -psi));

    // DT matrix
    arma::cx_mat DT(2, 2);
    DT(0, 0) = complex<double>(0.0, 0.0);
    DT(0, 1) = w;
    DT(1, 0) = complex<double>(0.0, 0.0);
    DT(1, 1) = complex<double>(0.0, 0.0);

    // DTC matrix
    arma::cx_mat DTC(2, 2);
    DTC(0, 0) = complex<double>(0.0, 0.0);
    DTC(0, 1) = complex<double>(0.0, 0.0);
    DTC(1, 0) = w;
    DTC(1, 1) = complex<double>(0.0, 0.0);

    // Print the matrices to verify
    LP.print("LP:");
    DT.print("DT:");
    DTC.print("DTC:");

    arma::cx_mat tunneling = arma::kron(w_on_g2, DTC) + arma::kron(w_on_g3, DTC) + arma::kron(arma::trans(w_on_g2), DT) + arma::kron(arma::trans(w_on_g3), DT) + arma::kron(arma::eye(37, 37), DT + DTC);
    arma::cx_mat layer_potential = arma::kron(Delta_Even, LP) + arma::kron(Delta_Odd, arma::conj(LP));
    arma::cx_mat potential = tunneling + layer_potential;
    // potential.print();

    // Initialize path_1_eig matrix
    mat path_1_eig(num_k, 74, fill::zeros); // Assuming 74 eigenvalues

    for (size_t i = 0; i < num_k; ++i)
    {
        double k_x = P1_x(i) / g_G;
        double k_y = P1_y(i) / g_G;
        cx_mat Kinetic = KE_SH(k_x, k_y, sixth_shell_mat);
        cx_mat Hamiltonian = potential + Kinetic;

        // Compute eigenvalues
        vec eigval = eig_sym(Hamiltonian);

        // Store the eigenvalues in the path_1_eig matrix
        path_1_eig.row(i) = eigval.t(); // Transpose to match the row assignment
    }

    mat path_2_eig(num_k, 74, fill::zeros); // Assuming 74 eigenvalues

    for (size_t i = 0; i < num_k; ++i)
    {
        double k_x = P2_x(i) / g_G;
        double k_y = P2_y(i) / g_G;
        cx_mat Kinetic = KE_SH(k_x, k_y, sixth_shell_mat);
        cx_mat Hamiltonian = potential + Kinetic;

        // Compute eigenvalues
        vec eigval = eig_sym(Hamiltonian);

        // Store the eigenvalues in the path_1_eig matrix
        path_2_eig.row(i) = eigval.t(); // Transpose to match the row assignment
    }

    mat path_3_eig(num_k, 74, fill::zeros); // Assuming 74 eigenvalues

    for (size_t i = 0; i < num_k; ++i)
    {
        double k_x = P3_x(i) / g_G;
        double k_y = P3_y(i) / g_G;
        cx_mat Kinetic = KE_SH(k_x, k_y, sixth_shell_mat);
        cx_mat Hamiltonian = potential + Kinetic;

        // Compute eigenvalues
        vec eigval = eig_sym(Hamiltonian);

        // Store the eigenvalues in the path_1_eig matrix
        path_3_eig.row(i) = eigval.t(); // Transpose to match the row assignment
    }

    mat path_4_eig(num_k, 74, fill::zeros); // Assuming 74 eigenvalues

    for (size_t i = 0; i < num_k; ++i)
    {
        double k_x = P4_x(i) / g_G;
        double k_y = P4_y(i) / g_G;
        cx_mat Kinetic = KE_SH(k_x, k_y, sixth_shell_mat);
        cx_mat Hamiltonian = potential + Kinetic;

        // Compute eigenvalues
        vec eigval = eig_sym(Hamiltonian);

        // Store the eigenvalues in the path_1_eig matrix
        path_4_eig.row(i) = eigval.t(); // Transpose to match the row assignment
    }

    path_4_eig.print("path_2_eig:");

    return 0;
}
