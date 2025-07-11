import numpy as np

def solve_y_realization(x_grid, epsilon, z):
    """
    Calculates y(x) for a given single realization of the random points z_i
    using the exact analytical solution.

    Args:
        x_grid (np.array): The grid of x-values to compute y on.
        epsilon (float): The small parameter epsilon.
        z (np.array): The sorted random locations of the delta functions.

    Returns:
        np.array: The solution y(x) on the specified grid.
    """
    L = 1.0 / epsilon
    e = np.exp(1.0)

    # Calculate the term needed for the derivative at the boundary y'(0)
    # Sum_i (e^(epsilon*(L-z_i)) - 1)
    sum_term_for_y_prime = np.sum(np.exp(epsilon * (L - z)) - 1)

    # Calculate y'(0)/epsilon from the boundary condition at x=L
    y_prime_0_over_eps = -(1.0 + epsilon * sum_term_for_y_prime) / (e - 1.0)

    # Calculate y(x) for each point in the grid
    y = np.zeros_like(x_grid)
    for i, x in enumerate(x_grid):
        # Find all z_j that are smaller than the current x
        z_less_than_x = z[z < x]
        
        # Calculate the sum over these points: Sum_{z_j<x} (e^(epsilon*(x-z_j)) - 1)
        if len(z_less_than_x) > 0:
            sum_term_for_y = np.sum(np.exp(epsilon * (x - z_less_than_x)) - 1)
        else:
            sum_term_for_y = 0.0

        # Combine the terms to get the final solution for y(x)
        term1 = 1.0
        term2 = y_prime_0_over_eps * (np.exp(epsilon * x) - 1)
        term3 = epsilon * sum_term_for_y
        
        y[i] = term1 + term2 + term3
        
    return y

def estimate_fluctuation_scaling(epsilons, num_realizations=500, num_x_points=200):
    """
    Estimates the fluctuation R for several epsilon values and finds the scaling law.
    """
    Rs = []
    print("Running simulations for each epsilon...")
    for epsilon in epsilons:
        L = 1.0 / epsilon
        N = int(L - 1)
        if N <= 0:
            print(f"Skipping epsilon={epsilon} as N <= 0.")
            continue
            
        x_grid = np.linspace(0, L, num_x_points)
        y_realizations = np.zeros((num_realizations, num_x_points))
        
        for i in range(num_realizations):
            # Generate N ordered uniform random values
            z = np.sort(np.random.uniform(0, L, N))
            y_realizations[i, :] = solve_y_realization(x_grid, epsilon, z)
            
        # The problem defines R based on Var[y(x) - y(0)].
        # Since y(0)=1 is a constant, Var[y(x)-y(0)] = Var[y(x)].
        var_y = np.var(y_realizations, axis=0)
        max_var = np.max(var_y)
        R = np.sqrt(max_var)
        Rs.append(R)
        print(f"  For epsilon = {epsilon:.3f}, R = {R:.5f}")
        
    # Fit a line to log(R) vs log(epsilon) to find the scaling
    log_eps = np.log(np.array(epsilons))
    log_Rs = np.log(np.array(Rs))
    
    # polyfit returns [slope, intercept]
    coeffs = np.polyfit(log_eps, log_Rs, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    
    print("\n--- Scaling Analysis ---")
    print("The relationship is assumed to be R = C * epsilon^p")
    print("Taking logs: log(R) = p * log(epsilon) + log(C)")
    print("\nFitting the simulation results to this line gives:")
    print(f"Slope (p): {slope:.4f}")
    print(f"Intercept (log(C)): {intercept:.4f}")
    
    C = np.exp(intercept)
    print(f"Pre-factor (C): {C:.4f}")
    
    print("\nFinal estimated scaling equation:")
    print(f"R(epsilon) ≈ {C:.4f} * epsilon^{slope:.4f}")

# Define the epsilon values to test
epsilons_to_test = [0.1, 0.05, 0.02, 0.01]
estimate_fluctuation_scaling(epsilons_to_test)
