import numpy as np

def estimate_limit(k_value, num_simulations):
    """
    Estimates the limit of n*P(n) for a given k via Monte Carlo simulation.

    Args:
        k_value (int): A large integer value for k. n is defined as 6*k.
        num_simulations (int): The number of simulations to run.
    """
    n = 6 * k_value
    # m is the number of vectors of each type, which is 2k.
    m = 2 * k_value
    
    # We can simulate the sum S by first simulating the sums of coefficients
    # for each vector type (N_A, N_B, N_C).
    # N_A is the sum of m Rademacher variables. This distribution is equivalent
    # to 2*B - m, where B is a binomial random variable B(m, 0.5).
    # This approach is much faster than simulating 6k variables for each trial.
    n_ones_A = np.random.binomial(m, 0.5, size=num_simulations)
    n_ones_B = np.random.binomial(m, 0.5, size=num_simulations)
    n_ones_C = np.random.binomial(m, 0.5, size=num_simulations)

    N_A = 2 * n_ones_A - m
    N_B = 2 * n_ones_B - m
    N_C = 2 * n_ones_C - m

    # The vectors are v_A=(1,0), v_B=(0.5, sqrt(3)/2), v_C=(-0.5, sqrt(3)/2).
    # The sum vector S = N_A*v_A + N_B*v_B + N_C*v_C.
    # Its components Sx and Sy are:
    Sx = N_A + 0.5 * N_B - 0.5 * N_C
    Sy = (np.sqrt(3.0) / 2.0) * (N_B + N_C)

    # The condition is ||S||_2 <= sqrt(2), which is ||S||_2^2 <= 2.
    norm_sq = Sx**2 + Sy**2
    
    # Count how many simulations satisfy the condition.
    successes = np.sum(norm_sq <= 2)

    # Estimate the probability P(n).
    P_n_estimated = successes / num_simulations

    # Estimate the value of n*P(n).
    n_Pn_estimated = n * P_n_estimated
    
    print(f"Running simulation with k = {k_value}, n = 6 * k = {n}")

    # As requested, printing numbers from the final analytical equation.
    # The final equation is derived from lim (6*k) / (3*k) as k -> infinity.
    analytic_numerator_factor = 6
    analytic_denominator_factor = 3
    analytic_limit = analytic_numerator_factor / analytic_denominator_factor

    print("\n--- Analytical Result ---")
    print("The limit can be derived using the Central Limit Theorem.")
    print(f"For large n=6k, P(n) is approximately 1/(3k).")
    print(f"The limit of n*P(n) is thus the limit of ({analytic_numerator_factor}*k)/({analytic_denominator_factor}*k).")
    print(f"Final equation: {analytic_numerator_factor} / {analytic_denominator_factor} = {analytic_limit}")

    print("\n--- Simulation Result ---")
    print(f"Number of simulations run: {num_simulations}")
    print(f"Number of sums S with ||S|| <= sqrt(2): {successes}")
    print(f"Estimated P(n) = {P_n_estimated:.8f}")
    print(f"Estimated n*P(n) = {n_Pn_estimated:.6f}")

if __name__ == '__main__':
    # Set parameters for the simulation.
    # A larger k and more simulations will yield a result closer to the analytical limit.
    k = 5000
    simulations = 2000000
    estimate_limit(k, simulations)
