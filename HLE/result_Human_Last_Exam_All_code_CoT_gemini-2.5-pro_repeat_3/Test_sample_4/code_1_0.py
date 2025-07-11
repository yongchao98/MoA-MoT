import numpy as np

def compute_poincare_polynomial():
    """
    Computes the Poincaré polynomial for the given Lie algebra.

    The calculation is based on the structural decomposition of the Lie algebra
    and properties of the Chevalley-Eilenberg cohomology.
    """

    # Coefficients of P_ker_d1(x) = 1 + x + x^2 + x^3, from lowest degree to highest.
    p_ker_d1 = [1, 1, 1, 1]

    # Coefficients of P_ker_d2(x) = 1 + x + x^2
    p_ker_d2 = [1, 1, 1]

    # Compute the product P_ker_d(x) = P_ker_d1(x) * P_ker_d2(x)
    # The result is the coefficients of the kernel's Poincaré polynomial.
    # numpy's polymul expects coefficients from highest degree to lowest, so we reverse.
    p_ker = np.polymul(p_ker_d1[::-1], p_ker_d2[::-1])[::-1]
    
    # For a degree-preserving derivation on a finite-dimensional graded vector space,
    # the Poincaré polynomials of the kernel and cokernel are identical.
    p_coker = p_ker

    # The Poincaré polynomial of g is P_g(x) = P_ker(x) + x * P_coker(x).
    # The coefficients of x * P_coker(x) are those of P_coker(x) shifted by one.
    p_coker_x = [0] + p_coker.tolist()

    # To add the polynomials, we make them the same length by padding with zeros.
    max_len = max(len(p_ker), len(p_coker_x))
    p_ker_padded = np.pad(p_ker, (0, max_len - len(p_ker)), 'constant')
    p_coker_x_padded = np.pad(p_coker_x, (0, max_len - len(p_coker_x)), 'constant')
    
    final_coeffs = p_ker_padded + p_coker_x_padded

    # Format the polynomial into a human-readable string.
    poly_str = []
    for i, coeff in enumerate(final_coeffs):
        if coeff == 0:
            continue
        
        # Handle the coefficient part
        if coeff == 1 and i > 0:
            coeff_str = ""
        else:
            coeff_str = str(int(coeff))

        # Handle the variable part
        if i == 0:
            var_str = ""
        elif i == 1:
            var_str = "x"
        else:
            var_str = f"x^{i}"
        
        # Add '*' for multiplication if needed
        if coeff_str and var_str:
            term = f"{coeff_str}*{var_str}"
        else:
            term = f"{coeff_str}{var_str}"

        poly_str.append(term)
        
    print(" + ".join(poly_str))

compute_poincare_polynomial()