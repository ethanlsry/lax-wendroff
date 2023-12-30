#include <iostream>
#include <cmath>
#include <complex>

std::vector<double> lax_wendroff(const std::vector<double> &u, double v, double dx, double dt){
    //use lax scheme to form intermediate u values
    std::vector<double> new_u = std::vector<double>(u.size());

    //iterate through vector except for ghost cells
    for (int j = 1; j < new_u.size()-1; j++){
        // std::cout << "j index: " << j << std::endl;
        double u_j_plus_1 = u[j+1];
        double u_j = u[j];
        double u_j_minus_1 = u[j-1];
        //use lax scheme to form intermediate values
        double u_n_plus_half_j_plus_half = 0.5 * (u_j + u_j_plus_1) - (((v * dt) / (2 * dx)) * (u_j_plus_1 - u_j));
        double u_n_plus_half_j_minus_half = 0.5 * (u_j_minus_1 + u_j) - ((v * dt) / (2 * dx)) * (u_j - u_j_minus_1);

        //use leapfrog scheme to compute u_j^{n+1}
        double u_j_n_plus_1 = u_j - ((v * dt) / (dx)) * (u_n_plus_half_j_plus_half - u_n_plus_half_j_minus_half);
        new_u[j] = u_j_n_plus_1;
    }
    
    //set ghost cells
    new_u[0] = u[u.size() - 2];
    new_u[u.size() - 1] = u[1];

    return new_u;
}