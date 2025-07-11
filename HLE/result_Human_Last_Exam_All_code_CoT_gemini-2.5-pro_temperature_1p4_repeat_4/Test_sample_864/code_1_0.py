import numpy as np

def solve_whitening_filter():
    """
    Solves for the whitening filter for a corrected, well-posed version of the problem.
    The corrected channel autocorrelation is Q(D) = 1*D^-1 + 2.5 + 1*D^1.
    """

    # Coefficients of the polynomial z*Q(z) = z^2 + 2.5*z + 1
    # The coefficients are ordered from highest power of z to the constant term.
    poly_coeffs = [1, 2.5, 1]

    # Find the roots of the polynomial
    roots = np.roots(poly_coeffs)
    # The roots are the zeros of Q(z)

    # In the D-transform convention where D=z^-1, a causal minimum-phase filter
    # has zeros |D| >= 1, which corresponds to zeros |z| <= 1.
    min_phase_zeros_z = [r for r in roots if np.abs(r) <= 1]
    
    # We will build the causal minimum-phase factor H(D).
    # Q(z) = (z - z1)(z - z2) where z1=-0.5, z2=-2.
    # Q(z) = (z+0.5)(z+2).
    # Q(D) = (1/D+0.5)(1/D+2) = (1+0.5D)/D * (1+2D)/D = (1+0.5D)(1+2D)/D^2.
    # H(D) should be causal (polynomial in D) and minimum-phase (zeros |D|>=1).
    # The zero at z=-0.5 corresponds to a zero at D = 1/z = -2. Since |-2|>=1, this is a min-phase zero.
    # The zero at z=-2 corresponds to a zero at D = 1/z = -0.5. Since |-0.5|<1, this is a max-phase zero.
    # So, the minimum-phase factor is related to (D - (-2)) = (D+2). Or more conventionally, (1 + D/2).
    # Let H(D) = C * (1 + 0.5*D)
    
    # We determine C by matching H(D)H(D^-1) to Q(D).
    # Q(D) = C^2 * (1 + 0.5*D) * (1 + 0.5*D^-1)
    # Q(D) = C^2 * (1 + 0.25 + 0.5*(D + D^-1))
    # Q(D) = C^2 * (1.25 + 0.5*(D + D^-1))
    # We want to match this to Q(D) = 2.5 + D + D^-1 = 2 * (1.25 + 0.5*(D + D^-1))
    # So, C^2 = 2 => C = sqrt(2).
    C = np.sqrt(2)
    a = 0.5 # The coefficient from (1 + a*D)

    # The causal, minimum-phase factor is H(D) = C * (1 + a*D)
    # H_D = C * (1 + a * D)

    # The anti-causal, maximum-phase factor is H(D^-1) = C * (1 + a*D^-1)
    
    # The whitening filter is W(D) = 1 / H(D^-1)
    # W(D) = 1 / (C * (1 + a*D^-1))
    
    c_w = 1 / C
    a_w = a

    print("The problem statement as given is ill-posed because the q_k sequence is not a valid autocorrelation sequence.")
    print("Solving for a corrected channel: Q(D) = 2.5 + D + D^-1.")
    print("\nThe causal, minimum-phase spectral factor of Q(D) is H(D) = C * (1 + a*D), where:")
    print(f"C = sqrt(2) = {C:.4f}")
    print(f"a = {a_w}")
    print("\nThe whitening filter W(D) is given by W(D) = 1 / H(D^-1).")
    print("This results in the following equation for W(D):")
    print("\nW(D) = c / (1 + a * D^-1)")
    print("\nWhere the calculated values for the coefficients are:")
    print(f"c = 1/C = {c_w:.4f}")
    print(f"a = {a_w}")
    print("\nFinal equation for the whitening filter:")
    print(f"W(D) = {c_w:.4f} / (1 + {a_w} * D^-1)")

solve_whitening_filter()