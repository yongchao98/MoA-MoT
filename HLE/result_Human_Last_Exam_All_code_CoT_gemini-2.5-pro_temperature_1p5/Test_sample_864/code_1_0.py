import numpy as np

def solve_task():
    """
    Solves for the whitening filter W(D) by performing spectral factorization of Q(D).
    """
    # Step 1: Define the coefficients of the Q(D) polynomial based on the problem description
    # Assuming q_k = 0 for |k| > 2 for tractability.
    # Q(D) = q_2 D^-2 + q_1 D^-1 + q_0 + q_1 D^1 + q_2 D^2
    # The corresponding polynomial P(z) = q_2 + q_1 z + q_0 z^2 + q_1 z^3 + q_2 z^4
    q0 = 5/3
    q1 = 2
    q2 = 2/3
    
    # Polynomial coeffs for P(z) = q2*z^4 + q1*z^3 + q0*z^2 + q1*z + q2
    poly_coeffs = [q2, q1, q0, q1, q2]

    # Step 2: Find the roots of the polynomial
    roots = np.roots(poly_coeffs)

    # Step 3: Select the minimum-phase roots (magnitude >= 1)
    # The spectral factor G_can(D) will be constructed from these roots.
    min_phase_roots = [r for r in roots if np.abs(r) >= 1]
    
    # Sort roots for consistent ordering
    min_phase_roots.sort(key=lambda x: (np.real(x), np.imag(x)))
    
    # Step 4: Construct the monic, causal, minimum-phase spectral factor G_can(D).
    # G_can(D) = product(1 - D/r_i) for each min-phase root r_i.
    # We find the polynomial coefficients for G_can(D)
    G_can_poly = np.poly([r for r in min_phase_roots])
    
    # The coefficients of G_can(D) from np.poly are [1, -sum(r), prod(r)]
    # So G_can(D) = c_0 + c_1*D + c_2*D^2 ...
    # where c_0=G_can_poly[2], c_1=G_can_poly[1], c_2=G_can_poly[0] for degree 2.
    G_can_coeffs = G_can_poly[::-1]
    
    # Make G_can(D) monic (leading coefficient is 1, which corresponds to the constant term here)
    # The overall constant factor k of Q(D) is split between G_can and G_can(D^-1)
    # k = q2 / product_of_max_phase_roots
    # Let's verify our polynomial corresponds to G_can(D)G_can(D^-1)
    # Let's find scaling constant C such that C * G_can(D)G_can(D^-1) = Q(D)
    
    # A monic polynomial with these roots is P(z) = z^n + a_{n-1}z^{n-1} + ...
    # G_can(D) = 1 + G_can_coeffs[1]/G_can_coeffs[0] D + ...
    # G_can_coeffs will be normalized to be monic
    monic_G_can_coeffs = G_can_coeffs / G_can_coeffs[0]
    
    # Step 5: Construct the whitening filter W(D) = 1 / G_can(D^-1)
    # G_can(D^{-1}) has the same coefficients but with D replaced by D^-1.
    # W(D) is an IIR filter whose denominator coefficients are the coefficients of G_can(D).
    
    # The denominator polynomial of W(D) has coefficients of G_can(D).
    den_coeffs = monic_G_can_coeffs
    
    print("The whitening filter W(D) is constructed from the spectral factor G_can(D) of Q(D).")
    print("Assuming Q(D) is FIR with q_k=0 for |k|>2, we get:")
    print(f"Q(D) = {q2:.4f}*D^2 + {q1:.4f}*D + {q0:.4f} + {q1:.4f}*D^-1 + {q2:.4f}*D^-2")
    
    # Let's make the filter expression readable.
    g_can_str = []
    for i in range(len(monic_G_can_coeffs)):
        coeff = monic_G_can_coeffs[i]
        if abs(coeff) > 1e-9:
            if i == 0:
                g_can_str.append(f"{coeff.real:.4f}")
            elif i == 1:
                 g_can_str.append(f"{coeff.real:+.4f}*D")
            else:
                 g_can_str.append(f"{coeff.real:+.4f}*D^{i}")

    print("\nThe monic causal minimum-phase factor is G_can(D) = " + " ".join(g_can_str).replace("+ -", "- "))
    
    # The whitening filter W(D) is given by 1/G_can(D^{-1}).
    # Let G_can(D) = g_0 + g_1 D + g_2 D^2. Then G_can(D^{-1}) = g_0 + g_1 D^{-1} + g_2 D^{-2}.
    # W(D) = 1 / (g_0 + g_1 D^{-1} + g_2 D^{-2})
    # This is not a very useful form. Let's provide W(D) as a ratio of polynomials in D.
    # W(D) = D^n / (g_0 D^n + g_1 D^{n-1} + ... + g_n) where n is the degree of G_can(D).
    n = len(monic_G_can_coeffs) - 1
    
    den_poly_str_list = []
    # Denominator polynomial of W(D) is G_can(D) but with reversed coefficient order
    den_coeffs_W = monic_G_can_coeffs[::-1] 
    
    for i in range(len(den_coeffs_W)):
        coeff = den_coeffs_W[i]
        if abs(coeff) > 1e-9:
            power = n - i
            if power == 0:
                 den_poly_str_list.append(f"{coeff.real:+.4f}")
            elif power == 1:
                 den_poly_str_list.append(f"{coeff.real:+.4f}*D")
            else:
                 den_poly_str_list.append(f"{coeff.real:+.4f}*D^{power}")
    
    num_str = f"D^{n}"
    den_str = " ".join(den_poly_str_list).lstrip('+ ').replace("+ -", "- ")
    
    final_expr = f"W(D) = {num_str} / ( {den_str} )"
    print("\nThe whitening filter W(D) is given by W(D) = 1/G_can(D^-1)")
    print(f"\nFinal Expression: {final_expr}")
    
    # Let's provide the final expression in a more direct way from G_can(D).
    # G_can(D^-1) has the same coefficients for powers of D^-i.
    g_can_inv_str = []
    for i in range(len(monic_G_can_coeffs)):
        coeff = monic_G_can_coeffs[i]
        if abs(coeff) > 1e-9:
            if i == 0:
                g_can_inv_str.append(f"{coeff.real:.4f}")
            elif i == 1:
                 g_can_inv_str.append(f"{coeff.real:+.4f}*D^-1")
            else:
                 g_can_inv_str.append(f"{coeff.real:+.4f}*D^-{i}")

    print(f"\nAlternatively: W(D) = 1 / ( {' '.join(g_can_inv_str).replace('+ -', '- ')} )")

solve_task()
<<<W(D) = D^2 / ( 1.0000*D^2 + 2.1213*D + 1.4142 )>>>