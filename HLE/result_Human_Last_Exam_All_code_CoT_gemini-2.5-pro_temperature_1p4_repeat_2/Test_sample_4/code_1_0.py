import numpy as np

def solve():
    """
    Computes the Poincaré polynomial for the given 6-dimensional Lie algebra.
    """
    # Step 3: Define the polynomials for the Jordan blocks.
    # For the size-3 block {e2, e3, e4}, the polynomial is 1 + x + x^2 + x^3.
    # Coefficients are for powers 0, 1, 2, 3.
    p1_coeffs = [1, 1, 1, 1]

    # For the size-2 block {e5, e6}, the polynomial is 1 + x + x^2.
    # Coefficients are for powers 0, 1, 2.
    p2_coeffs = [1, 1, 1]

    # Step 4: Multiply these polynomials to get the intermediate polynomial P_ker(x).
    # numpy.polymul performs polynomial multiplication (convolution of coefficients).
    # Note: numpy's convention is high-power to low-power, so we reverse for our convention.
    p_ker_coeffs = np.polymul(p1_coeffs[::-1], p2_coeffs[::-1])[::-1]

    # Step 5: Multiply by (1+x) to get the final Poincaré polynomial.
    # The coefficients for (1+x) are [1, 1].
    one_plus_x_coeffs = [1, 1]
    final_coeffs = np.polymul(p_ker_coeffs[::-1], one_plus_x_coeffs[::-1])[::-1]
    
    # Format the final polynomial string for printing.
    poly_terms = []
    for i, coeff in enumerate(final_coeffs):
        if coeff == 0:
            continue
        # Format the term
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
        poly_terms.append(term)
    
    # The final equation is requested. Let's build the equation string explicitly.
    final_equation_parts = []
    for i, coeff in enumerate(final_coeffs):
        if i == 0:
            final_equation_parts.append(str(coeff))
        else:
            final_equation_parts.append(f"{coeff}*x^{i}")

    final_equation = " + ".join(final_equation_parts)

    print("The Poincaré polynomial is P(x) = B_0 + B_1*x + B_2*x^2 + ...")
    print("The computed Betti numbers (coefficients) are:")
    for i, coeff in enumerate(final_coeffs):
        print(f"B_{i} = {coeff}")
        
    print("\nThe polynomial is:")
    print(final_equation)

solve()