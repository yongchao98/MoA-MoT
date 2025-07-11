import numpy as np

def calculate_trace_of_covariance():
    """
    This function calculates the trace of the covariance matrix for the given sampling procedure
    using an analytical approach.
    """
    # --- Problem Parameters ---
    d = 101
    alpha = 3.0
    beta = 2.0

    # The trace of the covariance matrix is given by the formula:
    # Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2

    # --- Step 1: Calculate E[||v||^2] ---
    # The vector v is obtained by applying a Householder reflection matrix R to d.
    # Reflections are orthogonal transformations, so they preserve vector norms.
    # Thus, ||v||^2 = ||d||^2.
    # The vector d is constructed such that its squared norm is:
    # ||d||^2 = ((a-b)/(a+b))^2 + (4ab/(a+b)^2) = (a^2-2ab+b^2+4ab)/(a+b)^2 = (a+b)^2/(a+b)^2 = 1.
    # Since ||v||^2 is always 1, its expectation is also 1.
    E_norm_v_sq = 1.0

    # --- Step 2: Calculate ||E[v]||^2 ---
    # Due to linearity of expectation, E[v] = R * E[d]. We first find E[d].
    
    # The first component d_1 = (a-b)/(a+b).
    # If a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta), then Y = a/(a+b) follows a Beta(alpha, beta) distribution.
    # We can write d_1 = a/(a+b) - b/(a+b) = Y - (1-Y) = 2*Y - 1.
    # The expectation of a Beta(alpha, beta) variable is alpha / (alpha + beta).
    E_Y = alpha / (alpha + beta)
    E_d1 = 2 * E_Y - 1

    # The other components d_i for i > 1 depend on c ~ N(0, I).
    # Due to the symmetry of the standard normal distribution, E[c_j/||c||] = 0.
    # Thus, E[d_i] = 0 for i > 1.
    # This means E[d] is a vector with E_d1 in the first position and zeros elsewhere.
    # E[d] = [E_d1, 0, ..., 0]^T = E_d1 * e_1.
    
    # Now we compute E[v] = R * E[d].
    # R is a reflection across the hyperplane orthogonal to u = v1 - v2.
    # v1 = e1 and v2 = 1_d, so u = [0, -1, ..., -1]^T.
    # The action of R on e1 is R*e1 = e1 - 2/(u^T*u) * u * (u^T*e1).
    # The dot product u^T * e1 is 0.
    # So, R*e1 = e1.
    # Therefore, E[v] = R * (E_d1 * e1) = E_d1 * (R * e1) = E_d1 * e1.
    # E[v] is the same as E[d].
    
    # The squared norm of E[v] is:
    # ||E[v]||^2 = ||E_d1 * e1||^2 = E_d1^2 * ||e1||^2 = E_d1^2.
    norm_E_v_sq = E_d1**2

    # --- Step 3: Final Calculation ---
    # Substitute the calculated values back into the formula.
    trace_cov = E_norm_v_sq - norm_E_v_sq
    
    # Print the final equation with the computed values
    print("The trace of the covariance matrix is calculated as E[||v||^2] - ||E[v]||^2.")
    print(f"The final equation is: {E_norm_v_sq} - {norm_E_v_sq} = {trace_cov}")

calculate_trace_of_covariance()