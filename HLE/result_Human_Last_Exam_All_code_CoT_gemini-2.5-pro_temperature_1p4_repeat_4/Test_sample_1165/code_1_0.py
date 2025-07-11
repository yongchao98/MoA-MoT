import numpy as np

def solve_for_y_fluctuation(epsilon, x_grid, z):
    """
    Calculates y(x) - y0(x) for a single realization of the random points z_i.
    This function computes the exact analytical solution for a given set of z_i.
    
    Args:
        epsilon (float): The parameter epsilon.
        x_grid (np.ndarray): The grid of x-points to evaluate the solution on.
        z (np.ndarray): The array of random points z_i.

    Returns:
        np.ndarray: The fluctuation y(x) - y0(x) on the x_grid.
    """
    L = 1.0 / epsilon
    
    # y(x) = 1 + C/ε * (exp(εx) - 1) + ε * Σ_{z_i<x} (exp(ε(x-z_i)) - 1)
    # where C = y'(0) is determined by the boundary condition y(L)=0.
    
    # Calculate y'(0), denoted as c0
    # exp(εL) = exp(1)
    sum_term = np.sum(np.exp(epsilon * (L - z)) - 1)
    numerator = 1 + epsilon * sum_term
    denominator = np.exp(epsilon * L) - 1
    c0 = - (epsilon / denominator) * numerator

    # Calculate y(x) for each x in the grid
    y_values = np.zeros_like(x_grid)
    for i, x in enumerate(x_grid):
        # Find z_j that are less than the current x
        z_less_than_x = z[z < x]
        
        # Calculate the sum part of the solution
        if len(z_less_than_x) > 0:
            sum_y = np.sum(np.exp(epsilon * (x - z_less_than_x)) - 1)
        else:
            sum_y = 0
            
        y_val = 1.0 + (c0 / epsilon) * (np.exp(epsilon * x) - 1) + epsilon * sum_y
        y_values[i] = y_val
        
    # The mean-field solution is y0(x) = 1 - epsilon * x
    y0_values = 1 - epsilon * x_grid
    
    # Return the fluctuation
    return y_values - y0_values

def run_simulation(epsilons, num_trials, num_x_points):
    """
    Runs a Monte Carlo simulation to estimate R(epsilon).
    
    Args:
        epsilons (list): A list of epsilon values to test.
        num_trials (int): The number of Monte Carlo trials for each epsilon.
        num_x_points (int): The number of points in the spatial grid.
        
    Returns:
        dict: A dictionary mapping each epsilon to its calculated R value.
    """
    R_results = {}
    
    print("Starting simulation...")
    for epsilon in epsilons:
        L = 1.0 / epsilon
        # N = ε⁻¹ - 1
        N = int(np.floor(L)) - 1
        if N <= 0:
            print(f"Skipping epsilon = {epsilon} because N <= 0.")
            continue
            
        x_grid = np.linspace(0, L, num_x_points)
        
        # Array to store fluctuation results from all trials
        y_f_trials = np.zeros((num_trials, num_x_points))
        
        for i in range(num_trials):
            # Generate N ordered uniform random variables on [0, L]
            z = np.sort(np.random.uniform(0, L, N))
            
            # Solve for the fluctuation and store it
            y_f = solve_for_y_fluctuation(epsilon, x_grid, z)
            y_f_trials[i, :] = y_f
            
        # Calculate the variance at each point x across the trials
        variance_y_f = np.var(y_f_trials, axis=0)
        
        # R^2 is the maximum of the absolute variance
        R_squared = np.max(np.abs(variance_y_f))
        
        # R is the square root
        R = np.sqrt(R_squared)
        R_results[epsilon] = R
        print(f"  - Epsilon = {epsilon:.4f}, N = {N}, R(epsilon) = {R:.6f}")
        
    return R_results

def main():
    # Parameters for the simulation
    epsilons = [0.1, 0.05, 0.02, 0.01, 0.005]
    num_trials = 2000
    num_x_points = 500

    # Run the simulation
    results = run_simulation(epsilons, num_trials, num_x_points)

    # Perform a log-log linear regression to find the scaling exponent
    valid_eps = sorted(results.keys())
    if len(valid_eps) < 2:
        print("\nNot enough data points for regression.")
        return
        
    log_eps = np.log(np.array(valid_eps))
    log_R = np.log(np.array([results[e] for e in valid_eps]))

    # slope is the scaling exponent, intercept is log of the prefactor
    slope, intercept = np.polyfit(log_eps, log_R, 1)
    prefactor = np.exp(intercept)

    print("\n--- Analysis ---")
    print("The relationship between R and epsilon is modeled as R(epsilon) = C * epsilon^k")
    print("A log-log plot of R vs. epsilon gives a line with slope k.")
    print("\n--- Result ---")
    print(f"The estimated function is: R(epsilon) ≈ {prefactor:.4f} * epsilon^{slope:.4f}")

if __name__ == '__main__':
    main()
