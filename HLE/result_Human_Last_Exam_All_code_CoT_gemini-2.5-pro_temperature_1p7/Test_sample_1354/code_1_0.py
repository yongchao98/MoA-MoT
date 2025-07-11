import numpy as np

def calculate_trace_of_covariance():
    """
    Calculates the trace of the covariance matrix based on the provided sampling procedure.
    
    The calculation is simplified by analytical insights to Tr(Cov(v)) = 1 - (E[(a-b)/(a+b)])^2.
    This function uses a Monte Carlo simulation to estimate the expectation term.
    """
    # Parameters from the problem description
    alpha = 3.0
    beta = 2.0
    theta = 1.0
    
    # Number of samples for the Monte Carlo simulation
    num_samples = 2000000

    # Step 1: Generate a large number of samples for a and b
    # a ~ gamma(alpha, theta)
    # b ~ gamma(beta, theta)
    # In numpy's gamma distribution, the second parameter is the scale, which corresponds to theta.
    a_samples = np.random.gamma(alpha, theta, num_samples)
    b_samples = np.random.gamma(beta, theta, num_samples)
    
    # Step 2: Calculate the term (a-b)/(a+b) for each sample
    # This corresponds to d_1 in the problem description
    d1_samples = (a_samples - b_samples) / (a_samples + b_samples)
    
    # Step 3: Estimate the expectation E[(a-b)/(a+b)] by taking the sample mean
    estimated_mean_d1 = np.mean(d1_samples)
    
    # Step 4: Calculate the squared norm of the expectation of d.
    # As derived in the plan, ||E[d]||^2 = (E[d_1])^2
    norm_sq_E_d = estimated_mean_d1**2
    
    # The term E[||v||^2] is known to be 1.
    E_norm_sq_v = 1.0

    # Step 5: Calculate the trace of the covariance matrix using the formula
    # Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2 = 1 - ||E[d]||^2
    trace_cov_v = E_norm_sq_v - norm_sq_E_d

    print("--- Monte Carlo Simulation Results ---")
    print(f"Number of samples used: {num_samples}")
    print(f"Estimated E[(a-b)/(a+b)]: {estimated_mean_d1}")
    print("\n--- Final Calculation ---")
    print("Trace(Cov(v)) = E[||v||^2] - ||E[v]||^2")
    print(f"              = {E_norm_sq_v} - ({estimated_mean_d1})^2")
    print(f"              = {E_norm_sq_v} - {norm_sq_E_d}")
    print(f"              = {trace_cov_v}")

calculate_trace_of_covariance()