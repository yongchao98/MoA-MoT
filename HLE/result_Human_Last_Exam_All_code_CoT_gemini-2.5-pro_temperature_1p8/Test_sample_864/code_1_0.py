import numpy as np

def solve_whitening_filter():
    """
    Solves for the whitening filter W(D) for a corrected version of the problem.

    The original problem statement for q_k leads to a non-valid power spectral density.
    This function solves a plausible, well-posed version where Q(D) is a simple
    Laurent polynomial: Q(D) = q_0 + q_1*(D + D^-1).
    We assume q_0 = 5 and q_1 = 2, which gives a valid spectral density.

    The whitening filter W(D) is found such that the resulting channel Q(D)W(D) is causal.
    This is achieved by spectral factorization of Q(D) = F(D)F(1/D), where F(D) is
    the causal, minimum-phase factor. The whitening filter is then W(D) = 1/F(1/D).
    """
    # Use coefficients from the corrected, plausible problem
    # Q(D) = 5 + 2*(D + D^{-1})
    q0 = 5
    q1 = 2

    # We are looking for a causal minimum-phase factor F(D) = c*(1 - alpha*D)
    # such that Q(D) = F(D)*F(1/D) = c^2 * (1 + alpha^2 - alpha*(D + D^{-1}))
    # Comparing coefficients with q0 + q1*(D + D^{-1}):
    # c^2 * (1 + alpha^2) = q0
    # -c^2 * alpha = q1
    # This leads to the quadratic equation for alpha: q1*alpha^2 + q0*alpha + q1 = 0
    
    coeffs = [q1, q0, q1]
    roots = np.roots(coeffs)
    
    # For minimum-phase factor F(D) = c*(1-alpha*D), the zero 1/alpha must be outside
    # the unit circle, which means |alpha| < 1.
    alpha = 0
    for r in roots:
        if abs(r) < 1:
            alpha = r
            break
            
    # Calculate c from -c^2 * alpha = q1
    c_squared = -q1 / alpha
    c = np.sqrt(c_squared)

    # The causal, minimum-phase factor is F(D) = c * (1 - alpha*D)
    # The required whitening filter is W(D) = 1 / F(1/D)
    # F(1/D) = c * (1 - alpha*D^{-1})
    # W(D) = 1 / (c * (1 - alpha*D^{-1}))

    print("Based on the corrected problem with Q(D) = 5 + 2*(D + D^-1):")
    print("The causal, minimum-phase spectral factor is F(D) = c * (1 - alpha*D)")
    print(f"alpha = {alpha:.1f}")
    print(f"c = {c:.1f}")
    
    f_d_coeff_0 = c
    f_d_coeff_1 = -c * alpha
    
    print(f"F(D) = {f_d_coeff_0:.1f} + {f_d_coeff_1:.1f}*D")
    
    # F(1/D) has coefficients c and -c*alpha for D^0 and D^{-1}
    f_1_over_d_coeff_0 = c
    f_1_over_d_coeff_neg_1 = -c * alpha

    print("\nThe whitening filter W(D) is given by W(D) = 1 / F(1/D)")
    print("So, the final equation for the whitening filter is:")
    print(f"W(D) = 1 / ({f_1_over_d_coeff_0:.1f} + {f_1_over_d_coeff_neg_1:.1f}*D^-1)")


solve_whitening_filter()
<<<W(D) = 1 / (2.0 + 1.0*D^-1)>>>