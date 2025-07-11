import numpy as np

def solve_pde():
    """
    Solves the given PDE numerically to find -u(0,1)/2.
    """
    # 1. Grid parameters
    L = 20.0
    dx = 0.1
    T = 1.0
    dt = 0.001
    Nx = int(2 * L / dx) + 1
    Nt = int(T / dt)
    x_grid = np.linspace(-L, L, Nx)
    center_index = Nx // 2

    # 2. Initial condition functions
    def u_initial(x):
        """u(x,0)"""
        # Add a small epsilon to avoid division by zero at x -> -inf for exp(x)
        return -2.0 + (1.0 - np.tanh(x)) / (np.exp(x) + 1.0)

    def ut_initial(x):
        """du/dt(x,0)"""
        tanh_x = np.tanh(x)
        sech_x = 1.0 / np.cosh(x)
        sech_half_x_sq = (1.0 / np.cosh(x / 2.0))**2
        return 0.25 * (tanh_x - 1.0) * sech_half_x_sq * (tanh_x - sech_x - 2.0)

    def P_nonlinear(u):
        """Nonlinear term (u-1)u(u+2)"""
        return (u - 1.0) * u * (u + 2.0)

    # 3. Initialization for the time-stepping scheme
    u_t0 = ut_initial(x_grid)
    u0 = u_initial(x_grid)

    # Numerically compute spatial derivatives at t=0
    ux0 = np.zeros_like(u0)
    uxx0 = np.zeros_like(u0)
    for j in range(1, Nx - 1):
        ux0[j] = (u0[j + 1] - u0[j - 1]) / (2.0 * dx)
        uxx0[j] = (u0[j + 1] - 2 * u0[j] + u0[j - 1]) / (dx**2)
    # Use one-sided differences for boundaries
    ux0[0] = (u0[1] - u0[0]) / dx
    ux0[-1] = (u0[-1] - u0[-2]) / dx
    uxx0[0] = (u0[2] - 2*u0[1] + u0[0]) / dx**2
    uxx0[-1] = (u0[-1] - 2*u0[-2] + u0[-3]) / dx**2


    # Compute u_tt(x,0) from the PDE
    utt0 = -8.0 * (u_t0 + u0 * ux0 - (1.0 / 8.0) * uxx0 - P_nonlinear(u0))

    # Initialize u_curr (at t=0) and u_prev (at t=-dt)
    u_curr = u0
    u_prev = u0 - dt * u_t0 + 0.5 * dt**2 * utt0
    
    # 4. Time-stepping loop
    for _ in range(Nt):
        u_next = np.zeros_like(u_curr)
        
        # Calculate for interior points
        for j in range(1, Nx - 1):
            ux = (u_curr[j + 1] - u_curr[j - 1]) / (2.0 * dx)
            uxx = (u_curr[j + 1] - 2 * u_curr[j] + u_curr[j - 1]) / (dx**2)
            
            spatial_and_nonlinear = u_curr[j] * ux - (1.0 / 8.0) * uxx - P_nonlinear(u_curr[j])
            
            # Rearranged finite difference formula:
            # u_next = A*u_prev + B*u_curr + C*spatial_and_nonlinear
            # where A = (1-4dt)/(1+4dt), B = 2/(1+4dt), C = 8dt^2/(1+4dt)
            # which comes from u_next*(1/(2dt)+1/(8dt^2)) = u_prev*(-1/(2dt)+1/(8dt^2)) + 2*u_curr/(8dt^2) + S
            
            coeff_prev = (1.0 / (8.0 * dt**2) - 1.0 / (2.0 * dt))
            coeff_curr = (2.0 / (8.0 * dt**2))
            coeff_next = (1.0 / (2.0 * dt) + 1.0 / (8.0 * dt**2))

            u_next[j] = (coeff_prev * u_prev[j] + coeff_curr * u_curr[j] + spatial_and_nonlinear) / coeff_next

        # Apply boundary conditions
        u_next[0] = 0.0
        u_next[-1] = -2.0
        
        # Update for the next iteration
        u_prev = u_curr.copy()
        u_curr = u_next.copy()

    # 5. Extract and print the final result
    u_at_0_1 = u_curr[center_index]
    final_result = -u_at_0_1 / 2.0
    
    print(f"The computed value of u(0, 1) is: {u_at_0_1}")
    print(f"The quantity -u(0,1)/2 is: -({u_at_0_1})/2 = {final_result}")

if __name__ == '__main__':
    solve_pde()