import numpy as np

def solve_trace_covariance():
    """
    Calculates the trace of the covariance matrix based on the analytical solution.
    """
    # Given parameters
    alpha = 3.0
    beta = 2.0
    
    # As explained in the plan, Tr(Cov(v)) simplifies to Tr(Cov(d)).
    # Tr(Cov(d)) = E[||d||^2] - ||E[d]||^2.

    # Step 1: Calculate E[||d||^2]
    # The vector d is constructed such that its norm ||d|| is always 1.
    # d_1 = (a-b)/(a+b)
    # d_tail = (2*sqrt(ab))/(||c||*(a+b)) * c
    # ||d||^2 = d_1^2 + ||d_tail||^2
    #         = ((a-b)/(a+b))^2 + (4*a*b / (||c||^2 * (a+b)^2)) * ||c||^2
    #         = (a^2 - 2ab + b^2 + 4ab) / (a+b)^2
    #         = (a+b)^2 / (a+b)^2 = 1.
    # Therefore, the expectation E[||d||^2] is 1.
    E_norm_d_sq = 1.0
    
    # Step 2: Calculate ||E[d]||^2
    # We need to find the expectation of the vector d, E[d].
    # E[d_tail] is a zero vector because E[c/||c||] is zero due to the symmetry of
    # the multivariate normal distribution N(0, Sigma).
    # So we only need to compute E[d_1].
    # E[d_1] = E[(a-b)/(a+b)] = E[2*a/(a+b) - 1] = 2*E[a/(a+b)] - 1.
    # For a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta), a/(a+b) ~ Beta(alpha, beta).
    # The mean of a Beta(alpha, beta) distribution is alpha / (alpha + beta).
    
    # E[a/(a+b)]
    E_beta = alpha / (alpha + beta)
    
    # E[d_1]
    E_d1 = 2 * E_beta - 1
    
    # E[d] = [E_d1, 0, 0, ..., 0].
    # So, ||E[d]||^2 = E_d1^2.
    norm_E_d_sq = E_d1**2
    
    # Step 3: Final Calculation
    # Tr(Cov(v)) = E[||d||^2] - ||E[d]||^2
    trace_cov = E_norm_d_sq - norm_E_d_sq

    print("--- Intermediate Calculations ---")
    print(f"Expectation of the first component of d, E[d_1]: {E_d1:.4f}")
    print(f"Squared norm of the mean vector d, ||E[d]||^2: {norm_E_d_sq:.4f}")
    print("\n--- Final Result ---")
    print(f"The trace of the covariance matrix is given by the equation: E[||d||^2] - ||E[d]||^2")
    print(f"Result: {E_norm_d_sq} - {norm_E_d_sq:.4f} = {trace_cov:.4f}")
    
    # Returning the exact value for the final answer block
    return trace_cov

# Execute the function to get the result
final_answer = solve_trace_covariance()
print(f"<<<{final_answer}>>>")
