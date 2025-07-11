import numpy as np

def solve_whitening_filter():
    """
    Calculates the whitening filter W(D) for the given ISI channel.
    The solution assumes the channel correlation sequence q_k is finite,
    terminating at k=4, based on the structure of the problem.
    """
    # Step 1: Define the coefficients q_k for k = 0, 1, 2, 3, 4.
    # q_k = 0 for |k| >= 5
    q = np.zeros(5)
    # k=0 (even)
    q[0] = 5/3 - 0/2
    # k=1 (odd)
    q[1] = 2 - (1-1)/2
    # k=2 (even)
    q[2] = 5/3 - 2/2
    # k=3 (odd)
    q[3] = 2 - (3-1)/2
    # k=4 (even)
    q[4] = 5/3 - 4/2
    
    # Laurent polynomial Q(D) has coefficients q_k for D^k
    # We construct a standard polynomial P(z) = z^4 * Q(z)
    # The coefficients of P(z) are [q_4, q_3, q_2, q_1, q_0, q_1, q_2, q_3, q_4]
    # from z^8 down to z^0
    poly_coeffs = np.array([q[4], q[3], q[2], q[1], q[0], q[1], q[2], q[3], q[4]])

    # Step 2: Find the roots of the polynomial P(z)
    roots = np.roots(poly_coeffs)

    # Step 3: Identify roots inside the unit circle for the minimum-phase factor G(D)
    min_phase_roots = roots[np.abs(roots) < 1]

    # Step 4: Construct the minimum-phase polynomial G(D)
    # G(D) = c * product(D - r_i) where r_i are the min_phase_roots
    # We can find c from q_4 = c^2 * product(r_i)
    prod_roots = np.prod(min_phase_roots)
    
    # We need to handle the sign of c. Let's construct the polynomial from roots
    # and then scale it.
    g_poly_unscaled = np.poly(min_phase_roots) # Corresponds to G(D)
    
    # G(D) = g_0 + g_1*D + ... + g_4*D^4
    # G(1/D) = g_0 + g_1/D + ... + g_4/D^4
    # Q(D) = G(D)*G(1/D).
    # The coefficient of D^4 in Q(D) is q_4.
    # The coefficient of D^4 in G(D)*G(1/D) is g_4*g_0
    # g_poly_unscaled is monic, so its g_4_unscaled is 1.
    # g_0_unscaled is (-1)^n * product(roots).
    g0_unscaled = g_poly_unscaled[-1]
    g4_unscaled = g_poly_unscaled[0]
    
    # q_4 = (c*g_4_unscaled) * (c*g_0_unscaled) = c^2 * g_4_unscaled * g_0_unscaled
    c_squared = q[4] / (g4_unscaled * g0_unscaled)
    c = np.sqrt(c_squared)
    
    # We must choose the sign of c. A common convention is to have g_0 > 0.
    if c * g_poly_unscaled[-1] < 0:
        c = -c
        
    g_poly = c * g_poly_unscaled
    
    # Step 5: The whitening filter is W(D) = 1/G(1/D)
    # G(1/D) = (g_0*D^4 + g_1*D^3 + g_2*D^2 + g_3*D + g_4) / D^4
    # So, W(D) = D^4 / (g_0*D^4 + g_1*D^3 + ... + g_4)
    
    # Denominator of W(D) has coefficients which are the reversed coefficients of G(D)
    w_den_coeffs = g_poly[::-1]
    # Numerator of W(D) is D^4
    w_num_coeffs = np.zeros(len(w_den_coeffs))
    w_num_coeffs[0] = 1.0

    print("The whitening filter W(D) is a rational function Num(D)/Den(D).")
    print("Where D is the unit delay operator.")
    print("\nNumerator Coefficients (from D^4 down to D^0):")
    for coeff in w_num_coeffs:
        print(f"{coeff:.4f}")
    
    print("\nDenominator Coefficients (from D^4 down to D^0):")
    for coeff in w_den_coeffs:
        print(f"{coeff:.4f}")
        
    # The final answer format is quite specific. It wants an equation.
    # W(D) = D^4 / (c0*D^4 + c1*D^3 + c2*D^2 + c3*D + c4)
    num_str = "D^4"
    den_str_parts = []
    for i, coeff in enumerate(w_den_coeffs):
        power = 4 - i
        if abs(coeff) > 1e-9: # a small tolerance for zero
            sign = "-" if coeff < 0 else "+"
            val = abs(coeff)
            if power > 1:
                den_str_parts.append(f"{sign} {val:.4f}*D^{power}")
            elif power == 1:
                den_str_parts.append(f"{sign} {val:.4f}*D")
            else:
                den_str_parts.append(f"{sign} {val:.4f}")

    # Clean up the start of the denominator string
    den_str = " ".join(den_str_parts).lstrip("+ ")

    final_equation = f"W(D) = {num_str} / ( {den_str} )"
    print("\nThe whitening filter W(D) is:")
    print(final_equation)

solve_whitening_filter()