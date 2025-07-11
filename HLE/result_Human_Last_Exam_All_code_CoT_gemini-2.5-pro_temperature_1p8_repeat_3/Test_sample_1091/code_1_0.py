import numpy as np

def solve_and_simulate():
    """
    This function provides the analytical solution to the problem and
    runs a Monte Carlo simulation to verify it.
    """

    # --- Analytical Solution ---
    print("--- Analytical Solution ---")
    
    # The limit of n * P(n) is given by 6 * f_Q(0)
    # where f_Q(q) is the probability density function of the limiting
    # quadratic form Q at q.
    
    # The eigenvalues of the matrix of the quadratic form are 0, 3/2, 3/2.
    # This leads to a distribution for Q that is a scaled Chi-squared
    # distribution with 2 degrees of freedom (an exponential distribution).
    # Its PDF at 0 is 1/3.
    
    c = 6  # from n = 6k
    f_Q_0 = 1/3 # PDF of Q at 0
    
    limit = c * f_Q_0
    
    print("The final equation for the limit is c * f_Q(0).")
    print(f"The constant c is: {c}")
    print(f"The value of the PDF at 0, f_Q(0), is: {f_Q_0:.4f}")
    print(f"The limit of n*P(n) as n -> infinity is {c} * {f_Q_0:.4f} = {limit:.4f}\n")


    # --- Numerical Simulation ---
    print("--- Monte Carlo Simulation ---")
    print("This simulation estimates n * P(n) for a large n to verify the analytical result.")

    # Parameters for the simulation
    # A larger k leads to a better approximation of the limit.
    k = 1000
    n = 6 * k
    num_trials = 100000
    
    two_k = 2 * k
    count_successful = 0

    print(f"Running simulation with k = {k} (n = {n}) and {num_trials} trials...")

    # We use vectorization for generating Rademacher variables to speed things up,
    # but process trials in a loop to keep memory usage low.
    for _ in range(num_trials):
        # Generate random signs (+1 or -1) for each sum
        # np.random.randint(0, 2, size=...) generates 0s and 1s
        # 2 * val - 1 maps {0, 1} to {-1, 1}
        eps_a = 2 * np.random.randint(0, 2, size=two_k) - 1
        eps_b = 2 * np.random.randint(0, 2, size=two_k) - 1
        eps_c = 2 * np.random.randint(0, 2, size=two_k) - 1
        
        A = np.sum(eps_a)
        B = np.sum(eps_b)
        C = np.sum(eps_c)
        
        # Calculate the squared norm of the sum S
        s_norm_sq = A**2 + B**2 + C**2 + A*B - A*C + B*C
        
        # Check if the norm is within the specified bound
        if s_norm_sq <= 2:
            count_successful += 1

    # Estimate P(n) = (number of successes) / (total trials)
    p_n_estimated = count_successful / num_trials
    
    # Calculate the final estimated value
    estimated_limit = n * p_n_estimated
    
    print(f"Simulation finished.")
    print(f"Number of successful trials (||S||_2^2 <= 2): {count_successful}")
    print(f"Estimated P(n) = {count_successful}/{num_trials} = {p_n_estimated:.6f}")
    print(f"Estimated value for n * P(n) is {n} * {p_n_estimated:.6f} = {estimated_limit:.4f}")
    
solve_and_simulate()