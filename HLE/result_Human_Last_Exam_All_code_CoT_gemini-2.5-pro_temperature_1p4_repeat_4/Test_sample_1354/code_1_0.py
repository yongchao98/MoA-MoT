import numpy as np

def solve_covariance_trace():
    """
    This function calculates the trace of the covariance matrix for the specified sampling procedure.
    It follows an analytical approach, explaining each step.
    """
    # Problem parameters
    d = 101
    alpha = 3.0
    beta = 2.0

    print("The trace of the covariance matrix is calculated using the formula:")
    print("Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2\n")

    # --- Part 1: Calculate E[||v||^2] ---
    print("Part 1: Calculating E[||v||^2]")
    print("The vector v is a Householder reflection of vector d, which is an orthogonal transformation.")
    print("Orthogonal transformations preserve the L2 norm, so ||v||^2 = ||d||^2.")
    print("The squared norm of d is ||d||^2 = ((a-b)/(a+b))^2 + (4ab)/(a+b)^2, which simplifies to 1.")
    E_norm_v_squared = 1.0
    print(f"Since ||v||^2 is always 1, its expectation E[||v||^2] = {E_norm_v_squared:.2f}\n")


    # --- Part 2: Calculate ||E[v]||^2 ---
    print("Part 2: Calculating ||E[v]||^2")
    # First, find E[d].
    # For a ~ Gamma(alpha,theta) and b ~ Gamma(beta,theta), we know that a/(a+b) ~ Beta(alpha, beta).
    # E[d_1] = E[(a-b)/(a+b)] = E[2*a/(a+b) - 1] = 2 * E[a/(a+b)] - 1
    E_a_div_apb = alpha / (alpha + beta)
    E_d1 = 2 * E_a_div_apb - 1
    print(f"The expectation of the first component of d is E[d_1] = 2 * ({alpha}/({alpha}+{beta})) - 1 = {E_d1:.2f}")

    # E[d_{2...d}] is 0 due to the symmetry of the distribution of c ~ N(0, I).
    print("The expectation of the other components of d is 0 due to symmetry.")
    E_d = np.zeros(d)
    E_d[0] = E_d1

    # Next, find E[v]. E[v] = H * E[d], where H is the Householder matrix.
    # The reflection vector is w = v1 - v2 = e1 - 1_d = [0, -1, ..., -1].
    # E[d] is proportional to e1. Since w is orthogonal to e1, w^T * E[d] = 0.
    # This means the reflection term is zero, so E[v] = E[d].
    E_v = E_d
    print("We find that E[v] = E[d] due to the orthogonality of the reflection vector and E[d].\n")

    # Now calculate ||E[v]||^2
    norm_E_v_squared = float(np.sum(E_v**2))
    print(f"Therefore, ||E[v]||^2 = ||E[d]||^2 = ({E_d1:.2f})^2 = {norm_E_v_squared:.2f}\n")


    # --- Part 3: Final Calculation ---
    print("Part 3: Final Calculation")
    trace_cov_v = E_norm_v_squared - norm_E_v_squared
    print("Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2")
    print(f"Plugging in the numbers: {E_norm_v_squared:.2f} - {norm_E_v_squared:.2f} = {trace_cov_v:.2f}")

solve_covariance_trace()