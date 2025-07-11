import numpy as np

def compute_poincare_polynomial():
    """
    Computes and prints the Poincaré polynomial for the given Lie algebra.
    """
    # The coefficients d_k = dim(ker(phi on wedge^k V))
    # where P_V(x) = sum(d_k * x^k)
    d = [1, 2, 4, 4, 2, 1]

    # The Poincaré polynomial of g is P_g(x) = (1+x) * P_V(x)
    # The coefficients of (1+x) are [1, 1]
    one_plus_x_coeffs = [1, 1]

    # Convolve the coefficient arrays to get the coefficients of the final polynomial
    betti_numbers = np.convolve(d, one_plus_x_coeffs)

    # Build the polynomial string representation
    poly_str = []
    for i, coeff in enumerate(betti_numbers):
        coeff = int(coeff)
        if coeff == 0:
            continue
        
        # Term part
        if i == 0:
            term = str(coeff)
        elif i == 1:
            if coeff == 1:
                term = "x"
            else:
                term = f"{coeff}*x"
        else:
            if coeff == 1:
                term = f"x^{i}"
            else:
                term = f"{coeff}*x^{i}"
        
        poly_str.append(term)

    print("The Poincaré polynomial is P(x) = " + " + ".join(poly_str))
    
    # Also printing the equation term by term as requested
    print("\nThe polynomial can be written as:")
    print("P(x) = " + "1 + 3*x + 6*x^2 + 8*x^3 + 6*x^4 + 3*x^5 + x^6")
    
compute_poincare_polynomial()