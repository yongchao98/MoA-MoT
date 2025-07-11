import numpy as np

def estimate_limit_nPn(k, num_trials):
    """
    Estimates the value of n*P(n) for a given k using Monte Carlo simulation.

    Args:
        k (int): Parameter determining the number of vectors. n = 6k.
        num_trials (int): The number of random sums to generate.
    """
    if k <= 0 or num_trials <= 0:
        print("k and num_trials must be positive integers.")
        return

    n = 6 * k
    
    # Define the three vectors
    v_A = np.array([1.0, 0.0])
    v_B = np.array([0.5, np.sqrt(3)/2])
    v_C = np.array([-0.5, np.sqrt(3)/2])
    
    # Counter for sums falling within the condition ||S||^2 <= 2
    count_in_disk = 0
    
    # Vectorized approach for efficiency
    # Generate all epsilons for all trials at once
    epsilons = 2 * np.random.randint(0, 2, size=(num_trials, n)) - 1
    
    # Sum epsilons for each vector type for each trial
    S_A_coeffs = np.sum(epsilons[:, 0 : 2*k], axis=1)
    S_B_coeffs = np.sum(epsilons[:, 2*k : 4*k], axis=1)
    S_C_coeffs = np.sum(epsilons[:, 4*k : 6*k], axis=1)
    
    # Calculate the sum vector S for each trial
    # S has shape (num_trials, 2)
    S_vectors = np.outer(S_A_coeffs, v_A) + np.outer(S_B_coeffs, v_B) + np.outer(S_C_coeffs, v_C)
    S_vectors = S_vectors.reshape(num_trials, 2)
    
    # Calculate the squared norm of each S vector
    squared_norms = np.sum(S_vectors**2, axis=1)
    
    # Count how many trials satisfy the condition
    count_in_disk = np.sum(squared_norms <= 2)

    # Estimate P(n)
    P_n = count_in_disk / num_trials
    
    # Calculate n * P(n)
    n_Pn = n * P_n
    
    print("--- Monte Carlo Simulation Results ---")
    print(f"Chosen k = {k}")
    print(f"Total number of vectors n = 6 * k = {n}")
    print(f"Number of trials = {num_trials}")
    print(f"Number of sums with ||S||^2 <= 2: {count_in_disk}")
    print(f"Estimated P(n) = {count_in_disk} / {num_trials} = {P_n:.6f}")
    print(f"Final estimated value of n*P(n) = {n} * {P_n:.6f} = {n_Pn:.4f}")
    print("The analytical result for the limit is 2.")

# --- Parameters for the simulation ---
# A larger k makes the CLT approximation more accurate.
# A larger num_trials makes the probability estimate more accurate.
k = 100
num_trials = 500000

# Run the simulation
estimate_limit_nPn(k, num_trials)
