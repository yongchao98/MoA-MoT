import numpy as np
from scipy.stats import linregress

def greens_function(x, s, L):
    """
    Green's function for y''=f(s) with y(0)=0, y(L)=0.
    """
    if x < s:
        return x * (s - L) / L
    else:
        return s * (x - L) / L

def estimate_R_for_epsilon(epsilon, num_runs=1000, num_x_points=200):
    """
    Estimates the value of R for a given epsilon using Monte Carlo simulation.
    """
    L = 1.0 / epsilon
    # N is the number of delta functions
    N = int(np.floor(1.0 / epsilon) - 1)
    if N <= 0:
        return 0

    x_grid = np.linspace(0, L, num_x_points)
    
    # Store the results of y2_stoch for each run
    y2_stoch_samples = np.zeros((num_runs, num_x_points))

    for i in range(num_runs):
        # Generate N random locations z_i, sorted as per problem description
        z_locations = np.sort(np.random.uniform(0, L, N))
        
        # Calculate y2_stoch(x) = sum over i of G(x, z_i)
        for j, x_val in enumerate(x_grid):
            g_sum = 0
            for z in z_locations:
                g_sum += greens_function(x_val, z, L)
            y2_stoch_samples[i, j] = g_sum
            
    # Calculate variance of y(x)
    # Var[y(x)] approx Var[epsilon^2 * y2_stoch(x)] = epsilon^4 * Var[y2_stoch(x)]
    var_y2 = np.var(y2_stoch_samples, axis=0)
    var_y = epsilon**4 * var_y2
    
    # Find the maximum variance and R
    max_var = np.max(var_y)
    R = np.sqrt(max_var)
    
    return R

def main():
    """
    Main function to run the simulation for several epsilon values and find the scaling.
    """
    # A set of epsilon values to test, on a log scale
    epsilon_values = np.logspace(-1, -2.5, 8)
    R_values = []

    print("Running simulations for various epsilon values...")
    for eps in epsilon_values:
        R = estimate_R_for_epsilon(eps)
        R_values.append(R)
        print(f"For epsilon = {eps:.4f}, estimated R = {R:.6f}")

    # Filter out any zero results
    valid_indices = [i for i, r in enumerate(R_values) if r > 0]
    if len(valid_indices) < 2:
        print("Not enough data to perform a linear regression.")
        return
        
    log_eps = np.log([epsilon_values[i] for i in valid_indices])
    log_R = np.log([R_values[i] for i in valid_indices])

    # Perform linear regression to find the scaling exponent
    slope, intercept, r_value, p_value, std_err = linregress(log_eps, log_R)
    
    C = np.exp(intercept)
    
    print("\n--- Scaling Analysis ---")
    print(f"Fit result: log(R) = {slope:.4f} * log(epsilon) + {intercept:.4f}")
    print("This corresponds to the scaling law: R(epsilon) = C * epsilon^k")
    print("The final equation is:")
    print(f"R(epsilon) = {C:.4f} * epsilon^{slope:.4f}")
    print(f"The theoretical scaling exponent is 0.5. The simulation gives {slope:.4f}.")

    print("\n--- Answer to the second question ---")
    print("Q: Do you expect the scaling for R(epsilon) to remain the same if z_i is an i.i.d. random variable, such that z_i ~ Normal(i, 0.5)?")
    print("A: Yes, the scaling R ~ epsilon^(1/2) is expected to remain the same.")
    print("The reasoning is that the fundamental components of the calculation do not change in terms of their epsilon dependence:")
    print("1. The number of sources N is still proportional to epsilon^(-1).")
    print("2. The sources z_i are still independent.")
    print("3. The Green's function, which represents the system's response to a single source, still has a magnitude that scales with L = epsilon^(-1).")
    print("The overall variance of the fluctuation is roughly N * (response_amplitude)^2. This leads to the same overall scaling with epsilon, although the constant pre-factor C would be different.")

if __name__ == '__main__':
    main()
