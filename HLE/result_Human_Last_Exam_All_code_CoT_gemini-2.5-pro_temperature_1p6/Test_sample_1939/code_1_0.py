import numpy as np
from scipy.special import binom

def solve_poly_properties(d):
    """
    Analyzes the Ehrhart polynomial p(z) = (z+1)^d for a given dimension d.
    
    Args:
    d (int): The dimension of the polytope's ambient space.
    """
    if d < 2:
        print("The problem is defined for d >= 2.")
        return

    # The Ehrhart polynomial is p(z) = (z+1)^d.
    # The coefficients of (z+1)^d are the binomial coefficients C(d, k).
    coeffs = [binom(d, k) for k in range(d + 1)]
    
    # In numpy, polynomials are represented with highest power first.
    poly = np.poly1d(coeffs)
    
    # Print the polynomial equation.
    # We construct the string representation manually for clarity.
    poly_str = "p(z) = "
    for i, coeff in enumerate(poly.c):
        power = d - i
        if coeff == 0:
            continue
        if i > 0:
            poly_str += " + "
        poly_str += f"{int(coeff)}*z^{power}"
    # Clean up the string representation
    poly_str = poly_str.replace("*z^1", "*z").replace("*z^0", "")
    print(f"For d={d}, the Ehrhart polynomial is {poly_str}")
    
    # Find the roots of the polynomial
    roots = poly.r
    print(f"\nThe roots of p(z) are: {roots}")

    # Evaluate the answer choices
    print("\n--- Checking the Answer Choices ---")
    
    # A. Every root of p has real part -1.
    is_A_correct = all(np.isclose(root.real, -1) for root in roots)
    print(f"A. Every root of p has real part -1.  -->  {is_A_correct}")

    # B. Every root of p is real.
    is_B_correct = all(np.isclose(root.imag, 0) for root in roots)
    print(f"B. Every root of p is real.  -->  {is_B_correct}")

    # C. The coefficients of p sum exactly d.
    coeff_sum = poly(1)
    is_C_correct = (coeff_sum == d)
    print(f"C. The coefficients of p sum to {coeff_sum}, not d={d}.  -->  {is_C_correct}")

    # D. The coefficients of p sum exactly d choose d/2.
    if d % 2 == 0:
        d_choose_d_half = binom(d, d // 2)
        is_D_correct = (coeff_sum == d_choose_d_half)
        print(f"D. The coefficients of p sum to {coeff_sum}, not C({d}, {d//2})={int(d_choose_d_half)}.  -->  {is_D_correct}")
    else:
        print(f"D. d/2 is not an integer for d={d}, so C(d,d/2) is not standardly defined for integers. --> False")


    # E. Every root of p has real part -1/2.
    is_E_correct = all(np.isclose(root.real, -0.5) for root in roots)
    print(f"E. Every root of p has real part -1/2.  -->  {is_E_correct}")
    
    print("\nConclusion: Statement A is the most specific correct answer.")

# We can run this for any dimension d >= 2, e.g., d=4
solve_poly_properties(d=4)