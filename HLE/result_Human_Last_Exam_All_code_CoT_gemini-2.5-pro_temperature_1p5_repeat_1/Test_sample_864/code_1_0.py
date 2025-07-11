import numpy as np

def solve_and_print_filter():
    """
    Solves for the whitening filter W(D) based on the problem description.
    """
    # Step 1: Define the polynomial P(D) from q_k
    # q_k for k = -4, -3, -2, -1, 0, 1, 2, 3, 4
    q_vals = np.array([
        -1/3., 1., 2/3., 2., 5/3., 2., 2/3., 1., -1/3.
    ])

    # P(D) = 3 * D^4 * Q(D)
    # The coefficients of Q(D) are q_k for D^k.
    # Q(D) = q[-4]D^-4 + ... + q[4]D^4
    # P(D) = 3 * (q[-4] + q[-3]D + ... + q[4]D^8)
    p_coeffs_rev = 3 * q_vals
    # numpy.roots expects coeffs from highest power to lowest
    p_coeffs = p_coeffs_rev[::-1]

    # Step 2: Find the roots of the polynomial P(D)
    roots = np.roots(p_coeffs)

    # Step 3: Identify roots for the minimum-phase factor G(D)
    # G(D) must have its roots outside or on the unit circle.
    # For roots on the unit circle, we take one from each conjugate pair.
    g_roots = []
    temp_roots = sorted(list(roots), key=lambda r: (np.abs(r), r.imag))
    
    processed_roots = set()
    for r in temp_roots:
        r_conj = np.conj(r)
        # Avoid processing a root if its conjugate has already been processed
        # (This handles the unit circle case where 1/r = conj(r))
        if any(np.isclose(r, pr) for pr in processed_roots):
            continue

        # For a reciprocal pair (r, 1/r*), pick the one outside the unit circle for G(D)
        if np.abs(r) > 1.000001:
            g_roots.append(r)
            processed_roots.add(r)
            processed_roots.add(1/r) # For real roots
            processed_roots.add(1/np.conj(r)) # For complex roots
        elif np.abs(r) > 0.999999: # Root is on the unit circle
            g_roots.append(r)
            processed_roots.add(r)
            processed_roots.add(r_conj)

    # Construct the monic polynomial G_monic(D) from its roots
    g_monic_coeffs = np.poly(g_roots)

    # Step 4: Determine the scaling constant C for G(D) = C * G_monic(D)
    # We use the relation q_4 = g_4 * g_0^*
    # g_4 = C (since G_monic is monic)
    # g_0 = C * G_monic(0) = C * (-1)^n * product_of_roots
    # q_4 = C * (C * (-1)^n * prod_roots)^* = |C|^2 * (-1)^n * conj(prod_roots)
    # where n is the degree of G(D) (which is 4)
    n = len(g_roots)
    prod_roots = np.prod(g_roots)
    q_n = q_vals[-1] # This is q_4 = -1/3
    
    C_squared = q_n / ((-1)**n * np.conj(prod_roots))

    if C_squared < 0:
        # This can happen due to numerical precision issues.
        # A small negative value should be treated as zero. Or it implies no real C exists.
        # The underlying theory guarantees Q(jw) >= 0, so C^2 must be real and >=0.
        # For simplicity, we take the real part and hope it's positive.
        C_squared = np.real(C_squared)
    
    C = np.sqrt(C_squared)
    g_coeffs = C * g_monic_coeffs

    # Step 5: The whitening filter W(D) is 1/G(1/D).
    # W(D) = D^n / G_rev(D) where G_rev has reversed coefficients of G.
    g_rev_coeffs = g_coeffs[::-1]

    # Print the result W(D) = Numerator(D) / Denominator(D)
    # Numerator is C' * D^4. Let's find the proper scaling.
    # W(D) = 1/G(1/D) = 1 / sum(g_k D^{-k}). We'll just present the denominator coefficients
    # Denominator polynomial is G(1/D) = g_0 + g_1 D^-1 + g_2 D^-2 + g_3 D^-3 + g_4 D^-4
    
    print("The problem as stated requires a numerical approach to factorize the channel correlation polynomial.")
    print("The causal minimum-phase factor G(D) of the channel correlation Q(D) is found to be:")
    
    g_expr = []
    for i, coeff in enumerate(g_coeffs):
        if not np.isclose(coeff, 0):
            g_expr.append(f"({coeff.real:.4f} + {coeff.imag:.4f}j)D^{n-i}")
    print(f"G(D) = " + " + ".join(g_expr).replace('+ -', '- '))
    
    print("\nThe required whitening filter W(D) that makes the channel causal is W(D) = 1 / G(1/D).")
    print("This results in the filter:")

    denom_poly_expr = []
    for i, coeff in enumerate(g_coeffs):
        if not np.isclose(coeff, 0):
            power = i
            term = f"({coeff.real:.4f} + {coeff.imag:.4f}j)"
            if power == 1:
                term += "D^{-1}"
            elif power > 1:
                term += f"D^{{-{power}}}"
            denom_poly_expr.append(term)
    
    print("W(D) = 1 / [" + " + ".join(denom_poly_expr).replace('+ -', '- ') + "]")

solve_and_print_filter()