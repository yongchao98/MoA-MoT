from fractions import Fraction

def solve_covariance_trace():
    """
    Analytically calculates the trace of the covariance matrix of the random variable v.
    The code follows the derivation explained in the plan.
    """

    # Parameters from the problem
    alpha = 3
    beta = 2

    # Step 1: Calculate E[d_1]
    # E[Z] for Z ~ Beta(alpha, beta) is alpha / (alpha + beta)
    # E[d_1] = 2 * E[Z] - 1
    # Using fractions for exact rational arithmetic
    frac_alpha = Fraction(alpha)
    frac_beta = Fraction(beta)
    
    frac_expected_Z = frac_alpha / (frac_alpha + frac_beta)
    frac_expected_d1 = 2 * frac_expected_Z - 1

    # Step 2: Calculate ||E[d]||^2
    # Since E[d] = [E[d_1], 0, ..., 0], the squared norm ||E[d]||^2 is just (E[d_1])^2
    frac_norm_sq_E_d = frac_expected_d1**2

    # Step 3: Calculate the final trace Tr(Cov(v)) = 1 - ||E[d]||^2
    frac_trace_cov_v = 1 - frac_norm_sq_E_d

    # Print the steps of the final calculation
    print("The trace of the covariance matrix, Tr(Cov(v)), simplifies to 1 - ||E[d]||^2.")
    print("\nFirst, we calculate the non-zero component of E[d]:")
    print(f"E[d_1] = 2 * (alpha / (alpha + beta)) - 1 = 2 * ({alpha}/({alpha}+{beta})) - 1 = {frac_expected_d1}")

    print("\nNext, we calculate the squared norm of E[d]:")
    print(f"||E[d]||^2 = (E[d_1])^2 = ({frac_expected_d1})^2 = {frac_norm_sq_E_d}")

    print("\nFinally, we compute the trace of the covariance matrix:")
    # Output the final equation with each number
    print(f"Tr(Cov(v)) = 1 - ||E[d]||^2 = 1 - {frac_norm_sq_E_d} = {frac_trace_cov_v}")
    
    print("\nThe numerical result is:")
    print(float(frac_trace_cov_v))

solve_covariance_trace()