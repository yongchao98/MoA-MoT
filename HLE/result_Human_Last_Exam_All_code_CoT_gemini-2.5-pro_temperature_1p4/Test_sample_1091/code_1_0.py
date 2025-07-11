import numpy as np

def solve_vector_sum_limit():
    """
    This function calculates the limit of n*P(n) as n->infinity based on the problem description.
    It prints the step-by-step derivation.
    """
    
    # We want to find the limit of n * P(n) as n -> infinity.
    # n = 6k
    # P(n) is the probability that the sum S of n random vectors satisfies ||S||_2 <= sqrt(2).

    # Step 1: Find the covariance matrix of the sum S.
    # S = (S_x, S_y). E[S] = 0.
    # The covariance matrix Sigma_S has components:
    # Var(S_x) = sum(v_ix^2) over all n=6k vectors
    # Var(S_y) = sum(v_iy^2) over all n=6k vectors
    # Cov(S_x, S_y) = sum(v_ix * v_iy) over all n=6k vectors
    # There are 2k vectors of type (1,0), 2k of type (0.5, sqrt(3)/2), and 2k of type (-0.5, sqrt(3)/2)
    # We can calculate the coefficients for k for each term in the covariance matrix.

    # Variance of the x-component of S, divided by k
    var_sx_k = 2 * (1.0)**2 + 2 * (0.5)**2 + 2 * (-0.5)**2
    # Variance of the y-component of S, divided by k
    var_sy_k = 2 * (0.0)**2 + 2 * (np.sqrt(3)/2)**2 + 2 * (np.sqrt(3)/2)**2
    # Covariance of S_x and S_y, divided by k
    cov_xy_k = 2 * (1.0)*(0.0) + 2 * (0.5)*(np.sqrt(3)/2) + 2 * (-0.5)*(np.sqrt(3)/2)

    # So, Sigma_S = [[var_sx_k * k, cov_xy_k * k], [cov_xy_k * k, var_sy_k * k]]
    
    # Step 2: Use the Central Limit Theorem to find P(n).
    # For large n, S follows a bivariate normal distribution with PDF f_S(x,y).
    # f_S(x,y) = (1 / (2*pi*sqrt(det(Sigma_S)))) * exp(-0.5 * r^T * Sigma_S^-1 * r)
    # At the center (0,0), f_S(0,0) = 1 / (2 * pi * sqrt(det(Sigma_S)))
    det_sigma_s_k2 = (var_sx_k * var_sy_k - cov_xy_k**2) # This is det(Sigma_S) / k^2
    sqrt_det_sigma_s_k = np.sqrt(det_sigma_s_k2) # This is sqrt(det(Sigma_S)) / k
    # f_S(0,0) = 1 / (2 * pi * sqrt_det_sigma_s_k * k)

    # P(n) is the integral of f_S over a disk with radius R = sqrt(2).
    # For large k, the PDF is almost constant over this disk, so P(n) ≈ Area * f_S(0,0).
    radius_squared = 2.0
    area = np.pi * radius_squared

    # P(n) ≈ area * f_S(0,0)
    # P(n) ≈ (pi * 2) * (1 / (2 * pi * 3 * k)) = 1 / (3k)

    # Step 3: Calculate the limit of n * P(n).
    # n = 6k, so k = n/6.
    # P(n) ≈ 1 / (3 * (n/6)) = 1 / (n/2) = 2/n
    # n * P(n) ≈ n * (2/n) = 2.

    # Now, let's print the derivation with the computed numbers.
    n_over_k = 6.0
    p_n_denominator_k_coeff = (2 * np.pi * sqrt_det_sigma_s_k) / area

    limit = n_over_k / p_n_denominator_k_coeff

    print("The calculation for the limit of n*P(n) as n -> infinity is as follows:")
    print("\n1. Compute the covariance matrix of the sum S:")
    print(f"Var(S_x) = 2k*(1)^2 + 2k*(0.5)^2 + 2k*(-0.5)^2 = {var_sx_k:.1f}k")
    print(f"Var(S_y) = 2k*(0)^2 + 2k*(sqrt(3)/2)^2 + 2k*(sqrt(3)/2)^2 = {var_sy_k:.1f}k")
    print(f"Cov(S_x, S_y) = 2k*(1*0) + 2k*(0.5*sqrt(3)/2) + 2k*(-0.5*sqrt(3)/2) = {cov_xy_k:.1f}k")
    print(f"So, the covariance matrix of S is Sigma_S = [[{var_sx_k:.1f}k, {cov_xy_k:.1f}k], [{cov_xy_k:.1f}k, {var_sy_k:.1f}k]].")

    print("\n2. Approximate the probability P(n) = P(||S||^2 <= 2):")
    print(f"The area of the disk where ||S|| <= sqrt(2) is pi * (sqrt(2))^2 = {radius_squared:.1f}*pi.")
    print(f"The PDF at S=0 is f(0) = 1 / (2*pi*sqrt(det(Sigma_S))) = 1 / (2*pi*{sqrt_det_sigma_s_k:.1f}k).")
    print(f"P(n) is approximated by Area * f(0):")
    print(f"P(n) approx = ({radius_squared:.1f}*pi) / (2*pi*{sqrt_det_sigma_s_k:.1f}k) = 1 / ({p_n_denominator_k_coeff:.1f}k).")
    
    print("\n3. Compute the final limit:")
    print("We are given n = 6k.")
    print(f"The expression for the limit is: lim (n -> inf) n * P(n)")
    print(f"= lim (k -> inf) ({n_over_k:.0f}k) * (1/({p_n_denominator_k_coeff:.0f}k))")
    print(f"= {n_over_k:.0f} / {p_n_denominator_k_coeff:.0f} = {limit:.0f}")

solve_vector_sum_limit()
<<<2>>>