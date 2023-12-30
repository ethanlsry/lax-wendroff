//use lax-wendroff function to solve 1D advection equation
#include "lax_wendroff.h"
#include <iostream>
#include <complex>
#include <fstream>

int main() {
    double x_min = 0.0;
    double x_max = 1.0;
    double n_x = 1000;
    double dx = (x_max - x_min) / n_x;
    double v = 1.0;
    double current_t = 0;
    double total_delta_t = 1.0;
    double dt = 0.9 * (dx / v);

    //calculate initial u
    std::vector<double> u_values_initial(n_x);
    double sigma = 0.05;
    double x = x_min;

    //set initial u values
    for (double i = x_min; i <= x_max+dx; i += dx) {
        int int_index = i*n_x;
        u_values_initial[int_index] = std::exp(-1 * ((((i - 0.5)*(i - 0.5)) / (sigma*sigma))));
    }

    //set ghost cells
    u_values_initial[0] = u_values_initial[n_x - 2];
    u_values_initial[n_x - 1] = u_values_initial[1];

    int iteration_tracker = 0;
    double print_at_this_time = 0.01;
    
    while (current_t <= total_delta_t+0.02){
        //calculate values at next delta t
        dt = 0.9 * (dx / v);

        //run it
        std::vector<double> u_values_next_time_step = lax_wendroff(u_values_initial, v, dx, dt);
        //update prior u values
        u_values_initial = u_values_next_time_step;

        //if current t is print time. print time is moving target incremented by 0.01
        if (current_t >= print_at_this_time){
            //save it
            std::string filename = "csv_output/data";
            filename += std::to_string(iteration_tracker);
            filename += ".csv";
            std::ofstream output_file(filename);
            for (double i = x_min; i <= x_max; i += dx) {
                int int_index = i*n_x;
                output_file << int_index << "," << u_values_next_time_step[int_index] << std::endl;
            }
            output_file.close();

            //increment target
            print_at_this_time += 0.01;
            iteration_tracker++;
        }
        // std::cout << "t=" << current_t << std::endl;

        //increment time
        current_t += dt;
    }
    return 0;
}