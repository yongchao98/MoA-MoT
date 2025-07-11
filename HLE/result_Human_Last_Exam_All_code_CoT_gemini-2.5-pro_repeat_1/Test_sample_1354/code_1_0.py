import numpy as np

def solve_trace_covariance():
    """
    Calculates the trace of the covariance matrix for the given sampling procedure.
    """
    # Parameters from the problem statement
    alpha = 3.0
    beta = 2.0
    d = 101

    # The trace of the covariance matrix is Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2.

    # Step 1: Calculate E[||v||^2]
    # The vector v is a reflection of the vector d. Reflections are orthogonal
    # transformations, so they preserve the norm. Thus, ||v|| = ||d||.
    # The vector d is constructed as d = [(a-b)/(a+b), 2*sqrt(ab)/(||c||*(a+b)) * c].
    # Its squared norm is ||d||^2 = ((a-b)/(a+b))^2 + (4ab/(a+b)^2) = (a^2-2ab+b^2+4ab)/(a+b)^2 = 1.
    # So, ||v||^2 is always 1.
    E_v_norm_sq = 1.0

    # Step 2: Calculate ||E[v]||^2
    # We have E[v] = H * E[d], where H is the reflection matrix.
    
    # First, we find E[d].
    # The first component is d_1 = (a-b)/(a+b).
    # Since a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta), the ratio Y = a/(a+b)
    # follows a Beta distribution, Y ~ Beta(alpha, beta).
    # E[Y] = alpha / (alpha + beta)
    # d_1 = a/(a+b) - b/(a+b) = Y - (1-Y) = 2*Y - 1.
    # E[d_1] = 2 * E[Y] - 1.
    E_Y = alpha / (alpha + beta)
    E_d1 = 2 * E_Y - 1

    # The other components of d involve the vector c ~ N(0, I_{d-1}).
    # Due to the spherical symmetry of the N(0,I) distribution, E[c/||c||] = 0.
    # Therefore, the expectation of the other components of d is 0.
    # So, E[d] is a vector with E_d1 in the first component and zeros elsewhere.
    # E[d] = E_d1 * e_1
    
    # Now, we compute E[v] = H * E[d] = H * (E_d1 * e_1) = E_d1 * (H * e_1).
    # The reflection matrix H is defined by u = v1 - v2.
    # v1 = e1 = [1, 0, ..., 0]^T
    # v2 = 1_d = [1, 1, ..., 1]^T
    # u = v1 - v2 = [0, -1, ..., -1]^T
    # The dot product u.T * e1 = 0, so e1 is orthogonal to u.
    # Vectors orthogonal to the reflection direction are unchanged by the reflection.
    # So, H * e1 = e1.
    # Therefore, E[v] = E_d1 * e1.
    
    # Now we can compute the squared norm of E[v].
    # ||E[v]||^2 = ||E_d1 * e1||^2 = E_d1^2 * ||e1||^2 = E_d1^2.
    norm_E_v_sq = E_d1**2

    # Step 3: Combine to find the trace
    trace = E_v_norm_sq - norm_E_v_sq

    # Output the final equation with all numbers
    print("The final equation is: Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2")
    print(f"E[||v||^2] = {E_v_norm_sq}")
    print(f"||E[v]||^2 = {norm_E_v_sq:.4f}")
    print(f"Tr(Cov(v)) = {E_v_norm_sq} - {norm_E_v_sq:.4f} = {trace:.4f}")

solve_trace_covariance()