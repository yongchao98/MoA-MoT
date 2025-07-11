import numpy as np
from scipy.optimize import curve_fit
import warnings

# Suppress harmless warnings for small epsilon
warnings.filterwarnings("ignore", category=RuntimeWarning)

def get_y_realization(x_grid, epsilon, z):
    """
    Calculates one realization of y(x) over a grid of x points using the analytical solution.
    """
    # The analytical solution for y(x) is derived by integrating the ODE.
    # It depends on y'(0), which itself depends on all the source locations z_i
    # to satisfy the boundary condition y(L)=0.
    
    # Calculate y'(0) for this realization
    # Note: L = 1/epsilon, so epsilon * L = 1
    sum_val = np.sum(np.exp(1.0 - epsilon * z) - 1.0)
    yp0 = -(epsilon / (np.e - 1.0)) * (1.0 + epsilon * sum_val)

    y_vals = np.zeros_like(x_grid)
    z_sorted = np.sort(z) # Sort z for efficient summation

    # Calculate y(x) for each x in the grid
    for i, x in enumerate(x_grid):
        # The term Σ_{z_i < x} is calculated efficiently by using the sorted z array
        relevant_z = z_sorted[z_sorted < x]
        
        sum_y = 0.0
        if len(relevant_z) > 0:
            sum_y = np.sum(np.exp(epsilon * (x - relevant_z)) - 1.0)
            
        y_vals[i] = 1.0 + yp0 * (np.exp(epsilon * x) - 1.0) / epsilon + epsilon * sum_y
        
    return y_vals

def estimate_R_for_epsilon(epsilon, n_sims=500, M=1000):
    """
    Estimates R for a given epsilon by running multiple simulations.
    """
    L = 1.0 / epsilon
    N = int(L - 1)
    if N <= 0:
        return np.nan

    x_grid = np.linspace(0, L, M)
    y_simulations = np.zeros((n_sims, M))

    for i in range(n_sims):
        # Generate N i.i.d. sources z_i uniformly distributed in (0, L)
        z_i_locations = np.random.uniform(0, L, N)
        y_simulations[i, :] = get_y_realization(x_grid, epsilon, z_i_locations)

    # Var[y(x) - y0(x)] = Var(y(x)) since y0 is deterministic.
    # We compute the variance of y(x) at each grid point across the simulations.
    variances = np.var(y_simulations, axis=0)
    
    # R is the square root of the maximum variance.
    max_variance = np.max(variances)
    R = np.sqrt(max_variance)
    return R

def power_law(x, C, alpha):
    """Power law function for fitting: R = C * epsilon^alpha"""
    return C * x**alpha

def main():
    """
    Main function to run the simulation, fit the results, and print the answer.
    """
    # Define a set of epsilon values to test
    epsilons = np.array([0.1, 0.05, 0.02, 0.01])
    
    # Calculate R for each epsilon
    print("Running simulations for different epsilon values...")
    R_values = [estimate_R_for_epsilon(eps) for eps in epsilons]
    print("Simulations finished.\n")
    
    # Fit the results to a power law
    try:
        params, _ = curve_fit(power_law, epsilons, R_values)
        C, alpha = params
        
        print("The estimated relationship between R and epsilon is of the form: R = C * epsilon^alpha")
        print("Based on the numerical simulation, the estimated parameters are:")
        print(f"C = {C:.4f}")
        print(f"alpha = {alpha:.4f}")
        
        print("\nThe final estimated equation is:")
        print(f"R(epsilon) = {C:.4f} * epsilon^{alpha:.4f}")

    except Exception as e:
        print(f"An error occurred during fitting: {e}")
        print("Could not determine the scaling relationship.")

if __name__ == "__main__":
    main()
