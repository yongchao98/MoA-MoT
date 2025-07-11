import numpy as np

def estimate_fluctuation_magnitude(epsilon, num_trials=2000, num_x_points=201):
    """
    Estimates the maximum magnitude of fluctuations for the given ODE.

    The function simulates the solution to the ODE for many random instances of
    the source term, calculates the variance, and returns the fluctuation
    magnitude R. It also prints the comparison with the analytical formula.
    """
    # Parameters from the problem description
    L = 1.0 / epsilon
    N = int(L) - 1

    # Grid for x variable
    x_grid = np.linspace(0, L, num_x_points)

    # To store the stochastic part of the solution for each trial
    # y(x) - y0(x) is approximated by epsilon^2 * S(x)
    # where S(x) = sum(G0(x, zi))
    S_values = np.zeros((num_trials, num_x_points))

    print(f"Starting simulation for epsilon = {epsilon}...")
    print(f"Domain L = {L}, Number of sources N = {N}")

    # Green's function for d^2/dx^2 with zero boundary conditions at 0 and L
    def green_function(x, z, L):
        if z <= x:
            return z * (L - x) / L
        else:
            return x * (L - z) / L

    # Run multiple trials to gather statistics
    for m in range(num_trials):
        # Generate N i.i.d. variables and sort them to get the order statistics z_i
        z = np.sort(np.random.uniform(0, L, N))

        # Calculate S(x_j) for this trial
        for j, x_j in enumerate(x_grid):
            g0_sum = sum(green_function(x_j, zi, L) for zi in z)
            S_values[m, j] = g0_sum

    # Calculate the variance of S(x) across trials for each x
    var_S = np.var(S_values, axis=0)

    # The variance of the fluctuation y(x) - y0(x) is epsilon^4 * Var(S)
    var_y_fluctuation = epsilon**4 * var_S

    # R is the sqrt of the maximum of this variance
    max_var = np.max(var_y_fluctuation)
    R_simulated = np.sqrt(max_var)

    # Analytical result for comparison
    R_analytical = np.sqrt(epsilon / 192.0)
    
    print("\n--- Results ---")
    print(f"Simulated R({epsilon}) = {R_simulated:.6f}")
    
    # Printing the components of the analytical equation as requested
    print("\nAnalytical Formula: R = sqrt(epsilon / 192)")
    print(f"sqrt({epsilon} / 192) = {R_analytical:.6f}")

# --- Main execution ---
# Set a value for epsilon
epsilon_val = 0.05
estimate_fluctuation_magnitude(epsilon_val)