import numpy as np

def calculate_trace_of_covariance():
    """
    This function calculates the trace of the covariance matrix for the given sampling procedure
    using a Monte Carlo simulation.
    """
    # --- Parameters ---
    d = 101
    alpha = 3.0
    beta = 2.0
    theta = 1.0  # Scale parameter for gamma distribution
    
    # As specified, v1 is the first standard basis vector and v2 is the vector of all ones.
    v1 = np.zeros(d)
    v1[0] = 1.0
    v2 = np.ones(d)
    
    # Number of samples for the Monte Carlo simulation
    N_samples = 200000

    # --- Pre-calculation for the Householder transformation ---
    # The transformation is v = (I - 2/(u.T*u) * u*u.T) @ d
    # which can be computed as v = d - (2 * (u.T @ d) / (u.T @ u)) * u
    u = v1 - v2
    # u.T @ u = ||u||^2. Since u = [0, -1, ..., -1], ||u||^2 = d - 1.
    u_norm_sq = float(d - 1)

    # --- Monte Carlo Simulation to estimate E[v] ---
    sum_v = np.zeros(d)
    for _ in range(N_samples):
        # Step 1: Sample a and b from gamma distributions
        a = np.random.gamma(alpha, theta)
        b = np.random.gamma(beta, theta)

        # Step 2: Sample c from a standard normal distribution N(0, I_{d-1})
        # Note: mu=0, Sigma=I_{d-1}
        c = np.random.randn(d - 1)

        # Step 3: Construct the vector d
        d_vec = np.zeros(d)
        a_plus_b = a + b
        
        # Guard against division by zero, although highly unlikely for continuous distributions
        if a_plus_b == 0:
            continue
            
        d_vec[0] = (a - b) / a_plus_b
        
        norm_c = np.linalg.norm(c)
        if norm_c == 0:
            # If c is the zero vector, the second component of d is undefined.
            # The probability of this is zero, so we can skip the sample.
            continue
        
        d_vec[1:] = (2.0 * np.sqrt(a * b) / (norm_c * a_plus_b)) * c

        # Step 4: Construct the vector v using the Householder reflection
        u_dot_d = np.dot(u, d_vec)
        v = d_vec - (2.0 / u_norm_sq) * u_dot_d * u
        
        # Accumulate the samples of v
        sum_v += v

    # --- Final Calculation ---
    # The trace is E[||v||^2] - ||E[v]||^2.
    # We showed analytically that E[||v||^2] = 1.
    E_norm_v_sq = 1.0
    
    # Estimate E[v] by averaging the samples
    E_v_hat = sum_v / N_samples
    
    # Estimate ||E[v]||^2
    norm_E_v_hat_sq = np.dot(E_v_hat, E_v_hat)
    
    # Calculate the final trace
    trace_cov_v = E_norm_v_sq - norm_E_v_hat_sq

    print("The trace of the covariance matrix is calculated as Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2.")
    print(f"From analytical derivation, E[||v||^2] = {E_norm_v_sq}")
    print(f"From simulation, the estimated ||E[v]||^2 is: {norm_E_v_hat_sq:.8f}")
    print(f"Final Result: {E_norm_v_sq} - {norm_E_v_hat_sq:.8f} = {trace_cov_v:.8f}")

calculate_trace_of_covariance()