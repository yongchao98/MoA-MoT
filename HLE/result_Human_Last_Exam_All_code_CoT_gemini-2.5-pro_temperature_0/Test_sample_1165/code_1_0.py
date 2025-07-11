import numpy as np
from scipy.stats import linregress

def estimate_fluctuation_scaling(epsilons, num_simulations, num_x_points=500):
    """
    Numerically estimates the scaling of fluctuation magnitude R with epsilon.

    Args:
        epsilons (list): A list of epsilon values to test.
        num_simulations (int): The number of simulations to run for each epsilon.
        num_x_points (int): The number of points to discretize the domain into.

    Returns:
        tuple: A tuple containing lists of epsilons and calculated R values.
    """
    calculated_Rs = []

    for epsilon in epsilons:
        L = 1.0 / epsilon
        N = int(L) - 1
        if N <= 0:
            print(f"Skipping epsilon={epsilon} as N={N} is not positive.")
            continue

        x_grid = np.linspace(0, L, num_x_points)
        
        # Accumulators for variance calculation: Var(X) = E[X^2] - (E[X])^2
        sum_y_minus_1 = np.zeros_like(x_grid)
        sum_y_minus_1_sq = np.zeros_like(x_grid)

        for _ in range(num_simulations):
            # Step 1: Generate ordered random z_i
            z = np.sort(np.random.uniform(0, L, N))

            # Step 2: Calculate the integration constant A for the exact solution
            # y(L) = 0
            # A = [-1 - ε * Σ(exp(1-εzᵢ)-1)] / [(e-1)/ε]
            sum_term = np.sum(np.exp(1 - epsilon * z) - 1)
            A_numerator = -1 - epsilon * sum_term
            A_denominator = (np.exp(1) - 1) / epsilon
            A = A_numerator / A_denominator

            # Step 3: Calculate y(x) across the grid for this simulation
            # y(x) = A/ε * (exp(εx)-1) + 1 + ε * Σ(exp(ε(x-zᵢ))-1)H(x-zᵢ)
            y_x = np.zeros_like(x_grid)
            for i, x in enumerate(x_grid):
                # The Heaviside sum is over z_i < x
                z_less_than_x = z[z < x]
                heaviside_sum = np.sum(np.exp(epsilon * (x - z_less_than_x)) - 1)
                
                y_x[i] = (A / epsilon) * (np.exp(epsilon * x) - 1) + 1 + epsilon * heaviside_sum
            
            # Step 4: Accumulate results for variance calculation
            y_minus_1 = y_x - 1.0
            sum_y_minus_1 += y_minus_1
            sum_y_minus_1_sq += y_minus_1**2

        # Step 5: Calculate variance and R for this epsilon
        mean_y_minus_1 = sum_y_minus_1 / num_simulations
        mean_y_minus_1_sq = sum_y_minus_1_sq / num_simulations
        variance = mean_y_minus_1_sq - mean_y_minus_1**2
        
        max_variance = np.max(variance)
        R = np.sqrt(max_variance)
        calculated_Rs.append(R)
        print(f"For ε = {epsilon:.3f}, R = {R:.5f}")

    return epsilons, calculated_Rs

if __name__ == '__main__':
    # Define simulation parameters
    # Use smaller epsilons to be in the asymptotic regime
    epsilons_to_test = [0.1, 0.05, 0.02, 0.01]
    simulations_per_epsilon = 1000

    # Run the simulation
    valid_epsilons, R_values = estimate_fluctuation_scaling(epsilons_to_test, simulations_per_epsilon)

    # Perform a linear fit on the log-log data to find the scaling exponent
    if len(valid_epsilons) > 1:
        log_eps = np.log(valid_epsilons)
        log_R = np.log(R_values)
        
        # The model is log(R) = p * log(ε) + log(C)
        slope, intercept, r_value, _, _ = linregress(log_eps, log_R)
        
        C = np.exp(intercept)
        p = slope

        print("\n--- Scaling Analysis ---")
        print(f"Fitting data to the model R(ε) = C * ε^p")
        print("The final equation with the computed numbers is:")
        print(f"R(ε) ≈ {C:.4f} * ε^{p:.4f}")
        print(f"The R-squared value of the fit is: {r_value**2:.4f}")
        print("The estimated scaling exponent p is the final answer.")
    else:
        print("Not enough data points to perform a scaling analysis.")
