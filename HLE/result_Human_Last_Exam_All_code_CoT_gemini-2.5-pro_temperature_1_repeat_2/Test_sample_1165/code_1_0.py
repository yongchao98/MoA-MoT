import numpy as np

def greens_function(x, z, L):
    """
    Calculates the Green's function G_0(x, z) for the operator d^2/dx^2
    on the interval [0, L] with zero boundary conditions.
    """
    if x <= z:
        return x * (z - L) / L
    else:
        return z * (x - L) / L

def estimate_R_and_print(epsilon, num_realizations=1000, num_x_points=201):
    """
    Performs a Monte Carlo simulation to estimate R(epsilon).
    
    Args:
        epsilon (float): The small parameter in the ODE.
        num_realizations (int): Number of Monte Carlo samples to generate.
        num_x_points (int): Number of points to discretize the domain [0, L].
    """
    if not (0 < epsilon < 1):
        print("Error: epsilon must be between 0 and 1.")
        return

    L = 1.0 / epsilon
    N = int(L) - 1

    if N <= 0:
        print(f"Epsilon ({epsilon}) is too large, resulting in N <= 0.")
        return

    print(f"Starting simulation for epsilon = {epsilon:.4f}")
    print(f"Domain L = {L:.1f}, Number of points N = {N}")
    print(f"Using {num_realizations} Monte Carlo realizations.")

    # Discretize the domain
    x_grid = np.linspace(0, L, num_x_points)
    
    # Store results for S(x) = sum_i G_0(x, z_i) for each realization
    s_values_all_realizations = np.zeros((num_realizations, num_x_points))

    # Run the Monte Carlo loop
    for k in range(num_realizations):
        # Generate N ordered random variables z_i from Uniform([0, L])
        z_i = np.sort(np.random.uniform(0, L, N))
        
        # Calculate S(x) for the current realization
        current_s_values = np.zeros(num_x_points)
        for i, x in enumerate(x_grid):
            s_val = sum(greens_function(x, z, L) for z in z_i)
            current_s_values[i] = s_val
        s_values_all_realizations[k, :] = current_s_values

    # Calculate the variance of S(x) at each point in the x_grid
    var_s = np.var(s_values_all_realizations, axis=0)

    # Find the maximum variance over the domain
    max_var_s = np.max(var_s)

    # Calculate R^2 = epsilon^4 * max_var_S
    R_squared = epsilon**4 * max_var_s
    R = np.sqrt(R_squared)

    # According to the theory, R(epsilon) = C * epsilon^(3/2)
    # Estimate the constant C
    C = R / (epsilon**1.5)

    print("\n--- Simulation Results ---")
    print(f"Maximum variance of the sum S(x) = {max_var_s:.4e}")
    print("The magnitude of fluctuations R is calculated as (max_x Var[y(x) - y0(x)])^0.5")
    
    # Print the final equation as requested
    print("\nFinal estimated equation:")
    print(f"R({epsilon:.4f}) = ({C:.4f}) * ({epsilon:.4f})^1.5 = {R:.4e}")
    print("--------------------------\n")


if __name__ == '__main__':
    # Set the value of epsilon for the simulation
    epsilon_value = 0.04  # This implies L=25, N=24
    estimate_R_and_print(epsilon_value)