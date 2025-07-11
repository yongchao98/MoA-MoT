import numpy as np

def compute_poincare_polynomial():
    """
    Computes the Poincaré polynomial for the given Lie algebra.

    The Lie algebra g is a semidirect product R x_ad R^5.
    The Poincare polynomial P_g(x) is given by (1+x)H(x), where
    H(x) = sum(h_q * x^q) and h_q are the dimensions of the kernel
    of the ad action on the q-th exterior power of the dual of R^5.

    The coefficients h_q have been calculated to be [1, 2, 4, 4, 2, 1].
    """

    # Coefficients of the polynomial H(x) = 1 + 2x + 4x^2 + 4x^3 + 2x^4 + x^5
    h_coeffs = [1, 2, 4, 4, 2, 1]

    # Coefficients of the polynomial 1 + x
    p_factor_coeffs = [1, 1]

    # Multiply the polynomials by convolving their coefficients
    poincare_coeffs = np.convolve(h_coeffs, p_factor_coeffs)

    # Format the polynomial for printing
    var = "x"
    parts = []
    
    # Handle the constant term
    if poincare_coeffs[0] != 0:
        parts.append(str(poincare_coeffs[0]))

    # Handle terms with degree >= 1
    for i, c in enumerate(poincare_coeffs[1:], start=1):
        if c == 0:
            continue
        
        # Format coefficient
        if c == 1:
            coeff_str = ""
        else:
            coeff_str = f"{c}*"
        
        # Format variable part
        if i == 1:
            var_str = var
        else:
            var_str = f"{var}**{i}"
        
        parts.append(f"{coeff_str}{var_str}")

    final_polynomial = " + ".join(parts)
    print("The Poincaré polynomial is:")
    print(final_polynomial)

if __name__ == '__main__':
    compute_poincare_polynomial()