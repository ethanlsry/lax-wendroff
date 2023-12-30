<h3>Solving PDEs via Lax-Wendroff</h3>

The Lax-Wendroff method solves a partial differential equation by advancing through timesteps and evaluating the system at each timestep. At each step the system is evaluated at half tiem steps and half grid points, then at a second step the system is evaluated at using prior time step values. This method is second order both spatially and temporally. Compared to other PDE numerical techniques such as the Finite Difference method, the Lax-Wendroff method can avoid erronious oscillations.

Created with support from Alex Chen.