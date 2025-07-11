from fractions import Fraction

def solve_trace_covariance():
    """
    This function calculates the trace of the covariance matrix based on the analytical derivation.
    """
    # Given parameters
    d = 101
    alpha = 3
    beta = 2
    theta = 1
    # v1 = e1, v2 = 1_d

    print("This script calculates the trace of the covariance matrix for the specified random variable v.")
    print("The analytical solution follows the formula: Tr(Cov(v)) = 1 - ||E[v]||^2\n")

    # Step 1: Calculate E[d_1]
    print("Step 1: Calculate the expectation of the first component of d, E[d_1].")
    # Using the property that a/(a+b) for a~Gamma(alpha, th), b~Gamma(beta, th) is Beta(alpha, beta)
    # The mean of Beta(alpha, beta) is alpha / (alpha + beta)
    E_X_frac = Fraction(alpha, alpha + beta)
    E_X_float = float(E_X_frac)
    
    print(f"Let X = a/(a+b). Given a ~ Gamma({alpha}, {theta}) and b ~ Gamma({beta}, {theta}), X follows a Beta({alpha}, {beta}) distribution.")
    print(f"E[X] = alpha / (alpha + beta) = {alpha} / ({alpha} + {beta}) = {E_X_frac}")

    # E[d_1] = E[2*X - 1] = 2*E[X] - 1
    E_d1_frac = 2 * E_X_frac - 1
    E_d1_float = float(E_d1_frac)
    print(f"d_1 = (a-b)/(a+b) = 2*X - 1. So, E[d_1] = 2 * E[X] - 1 = {E_d1_frac} (or {E_d1_float:.4f}).")
    print("-" * 20)

    # Step 2: Determine E[d]
    print("Step 2: Determine the expectation of the vector d, E[d].")
    print("For components i > 1, E[d_i] = 0 due to the symmetry of the N(0, I) distribution of c.")
    print(f"Therefore, E[d] = [{E_d1_frac}, 0, ..., 0]^T.")
    print("-" * 20)
    
    # Step 3: Determine E[v]
    print("Step 3: Determine the expectation of the vector v, E[v].")
    print("E[v] = H * E[d], where H is the Householder matrix.")
    print("The vector u = v1 - v2 gives u^T * e1 = 0, which simplifies H * e1 to e1.")
    print(f"Thus, E[v] = E[d] = [{E_d1_frac}, 0, ..., 0]^T.")
    print("-" * 20)

    # Step 4: Calculate ||E[v]||^2
    print("Step 4: Calculate the squared norm of E[v].")
    norm_E_v_sq_frac = E_d1_frac**2
    norm_E_v_sq_float = float(norm_E_v_sq_frac)
    print(f"||E[v]||^2 = ||E[d]||^2 = (E[d_1])^2 = ({E_d1_frac})^2 = {norm_E_v_sq_frac} (or {norm_E_v_sq_float:.4f}).")
    print("-" * 20)

    # Step 5: Calculate the final trace
    print("Step 5: Calculate the final trace Tr(Cov(v)).")
    trace_cov_v_frac = 1 - norm_E_v_sq_frac
    trace_cov_v_float = float(trace_cov_v_frac)
    print(f"The final equation is: Tr(Cov(v)) = 1 - ||E[v]||^2")
    print(f"Tr(Cov(v)) = {1} - {norm_E_v_sq_frac} = {trace_cov_v_frac}")
    print(f"The numerical result is: {trace_cov_v_float:.4f}")

    return trace_cov_v_float

if __name__ == '__main__':
    solve_trace_covariance()
