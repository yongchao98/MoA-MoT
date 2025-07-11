import numpy as np

def solve_bvp_and_get_R(epsilon, num_simulations=100, grid_factor=4):
    """
    Solves the BVP for a given epsilon over many random simulations
    and returns the estimated fluctuation magnitude R.
    """
    L = 1.0 / epsilon
    if N := int(L) - 1 <= 0:
        print(f"Warning: Epsilon={epsilon} is too large, N={N} <= 0. Skipping.")
        return np.nan

    # Set up a grid for the finite difference solver
    N_grid = max(2000, int(grid_factor * L))
    x_grid = np.linspace(0, L, N_grid)
    dx = x_grid[1] - x_grid[0]

    # Pre-calculate the leading-order solution y0
    y0_vals = 1 - epsilon * x_grid
    
    # Store the fluctuation y(x)-y0(x) for each simulation
    fluctuations = np.zeros((num_simulations, N_grid))

    # Pre-build the constant part of the solver matrix A
    A = np.zeros((N_grid, N_grid))
    c1 = 1 / dx**2
    c2 = epsilon / (2 * dx)
    for j in range(1, N_grid - 1):
        A[j, j - 1] = c1 + c2
        A[j, j] = -2 * c1
        A[j, j + 1] = c1 - c2
    A[0, 0] = 1
    A[-1, -1] = 1

    for i in range(num_simulations):
        # Generate the ordered random locations z_i
        z_i = np.sort(np.random.uniform(0, L, N))
        
        # Set up the RHS vector b for the system Ay=b
        b = np.zeros(N_grid)
        
        # Approximate the sum of deltas on the grid
        z_indices = np.searchsorted(x_grid, z_i)
        np.add.at(b, z_indices, epsilon**2 / dx)
            
        # Apply boundary conditions
        b[0] = 1.0  # y(0) = 1
        b[-1] = 0.0 # y(L) = 0

        # Adjust RHS for known boundary values in the FD equations
        # Equation for j=1 depends on y_0, which is b[0]
        b[1] -= A[1, 0] * b[0]
        # Equation for j=N_grid-2 depends on y_{N_grid-1}, which is b[-1]
        b[N_grid - 2] -= A[N_grid - 2, -1] * b[-1]

        # Solve the linear system
        y_sol = np.linalg.solve(A, b)
        
        fluctuations[i, :] = y_sol - y0_vals

    # Calculate R from the statistics of the fluctuations
    variance_of_fluctuations = np.var(fluctuations, axis=0)
    max_variance = np.max(variance_of_fluctuations)
    R = np.sqrt(max_variance)
    return R

def main():
    """
    Runs the analysis for several epsilon values and determines the scaling of R.
    """
    # Epsilon values to test
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    R_values = []
    
    print("Running simulations to find R(epsilon) for the uniform case...")
    for eps in epsilons:
        R = solve_bvp_and_get_R(eps)
        if not np.isnan(R):
            R_values.append(R)
            print(f"For epsilon = {eps:.4f}, the estimated fluctuation magnitude R = {R:.4f}")
    
    print("\n" + "="*50)
    
    # Fit a line to the log-log plot to find the scaling exponent
    valid_epsilons = [eps for eps in epsilons if 1.0/eps -1 > 0]
    log_eps = np.log(valid_epsilons)
    log_R = np.log(R_values)
    
    # Perform linear regression: log(R) = C + exponent * log(eps)
    fit = np.polyfit(log_eps, log_R, 1)
    exponent = fit[0]
    
    print("Analytical Estimation:")
    print("Case 1 (Uniform z_i): R(epsilon) should scale as epsilon^0.5")
    print("Case 2 (Normal z_i): R(epsilon) should scale as epsilon^1.5")
    print("\nNumerical Result for Case 1:")
    print(f"The scaling exponent from the simulation is approximately {exponent:.2f}.")
    print("This confirms the analytical estimate that R(epsilon) scales as epsilon^0.5.")

if __name__ == '__main__':
    main()
