import numpy as np

def solve_whitening_filter():
    """
    Solves for the whitening filter W(D) for the corrected ISI channel.
    """
    # Original coefficients q_k for k = 0, 1, 2, 3, 4
    q = np.array([5./3., 2., 2./3., 1., -1./3.])
    
    # Laurent polynomial coefficients for Q(D) from D^-4 to D^4
    Q_coeffs = np.array([q[4], q[3], q[2], q[1], q[0], q[1], q[2], q[3], q[4]])
    
    # At D = -1, Q(-1) = q0 - 2*q1 + 2*q2 - 2*q3 + 2*q4
    q_neg1_val = q[0] - 2*q[1] + 2*q[2] - 2*q[3] + 2*q[4]
    
    print("The D-transform Q(D) is constructed from the given q_k.")
    print(f"Evaluating at D=-1, we get Q(-1) = {q_neg1_val:.4f}, which is negative.")
    print("This means Q(D) is not a valid power spectrum, so we correct it by adding a constant.")
    
    # Corrected Q'(D) = Q(D) - Q(-1). This is equivalent to changing q0.
    Q_prime_coeffs = np.copy(Q_coeffs)
    Q_prime_coeffs[4] -= q_neg1_val
    
    print(f"The corrected spectrum is Q'(D) = Q(D) - ({q_neg1_val:.4f}) = Q(D) + {abs(q_neg1_val):.4f}.")
    
    # Convert to a standard polynomial P(z) = z^4 * Q'(z) to find roots
    # P(z) has coefficients from z^0 to z^8
    P_coeffs = np.flip(Q_prime_coeffs) 
    
    # Find roots of the polynomial P(z)
    roots = np.roots(P_coeffs)
    
    # Find the causal, minimum-phase factor H_c(D)
    # Its roots must be outside the unit circle.
    # We construct H_c(D) = C * product(D - r_i) for roots |r_i| > 1
    # For a monic H_c(D) (h_0 = 1), we use H_c(D) = product(1 - D/r_i)
    
    roots_outside = roots[np.abs(roots) > 1.00001]
    roots_on = roots[np.abs(np.abs(roots) - 1.0) < 0.00001]
    
    # Hc(D) is built from roots outside and on the unit circle
    # Let's collect all roots for our Hc(D). For pairs (r, 1/r), we take the one outside.
    # For pairs on the unit circle (r, 1/r=conj(r)), we take one of them.
    selected_roots = []
    temp_roots = list(roots)
    while len(temp_roots) > 0:
        r1 = temp_roots.pop(0)
        # Find reciprocal root to remove it
        found_reciprocal = False
        for i, r2 in enumerate(temp_roots):
            if np.isclose(r1, 1/r2):
                temp_roots.pop(i)
                found_reciprocal = True
                break
        
        if np.abs(r1) > 1:
            selected_roots.append(r1)
        elif found_reciprocal: # r1 is inside, its reciprocal is outside
            selected_roots.append(1/r1)
        else: # root is on the unit circle
            selected_roots.append(r1)

    # Construct the monic minimum phase Hc(D) = product(1 - D/r_i)
    # The coefficients of this polynomial are found by convolving (1, -1/r_i) terms
    hc_coeffs = np.array([1.])
    for r in selected_roots:
        hc_coeffs = np.convolve(hc_coeffs, np.array([1., -1./r]))
    
    print("\nThe spectral factor H_c(D) is found. It is causal and minimum-phase.")
    # The whitening filter W(D) = 1 / H_c(1/D)
    # H_c(1/D) has coefficients of H_c(D) in reverse order.
    # So W(D) is a rational function.
    # Denominator of W(D) is poly with coefficients from hc_coeffs in reverse order.
    # Numerator of W(D) is 1.
    
    w_num_coeffs = np.array([1.0])
    w_den_coeffs = np.flip(hc_coeffs) # These are coeffs of Hc(1/D)

    print("\nThe whitening filter W(D) is given by a ratio of polynomials, W(D) = N(D)/M(D).")
    print("Numerator Coefficients N(D):")
    for c in w_num_coeffs:
        print(f"{c.real:.4f} ", end="")
    print("\nDenominator Coefficients M(D) (for powers D^0, D^1, ...):")
    for c in w_den_coeffs:
        print(f"{c.real:.4f} ", end="")
    print("\n")
    
    # The final equation is Q_eq(D) = H_c(D)
    print("The equalized channel response Q_eq(D) = H_c(D) has the following coefficients (for D^0, D^1, ...):")
    final_eq = ""
    for i, c in enumerate(hc_coeffs):
        if i == 0:
             final_eq += f"{c.real:.4f}"
        else:
             sign = '+' if c.real >= 0 else '-'
             final_eq += f" {sign} {abs(c.real):.4f}*D^{i}"
    print(final_eq)

solve_whitening_filter()