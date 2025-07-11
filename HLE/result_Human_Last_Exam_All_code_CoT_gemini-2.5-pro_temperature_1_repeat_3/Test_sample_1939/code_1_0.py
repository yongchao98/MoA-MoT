import numpy as np

def solve_poly():
    """
    This script analyzes the Ehrhart polynomial for the specified polytope for the case d=3.
    The polynomial has been determined to be p(z) = 0.5*z^3 + 1.5*z^2 + 2*z + 1.
    We will find the roots of this polynomial to check the given answer choices.
    To work with integer coefficients, we analyze 2*p(z) = z^3 + 3*z^2 + 4*z + 2.
    """
    
    # Coefficients of the polynomial z^3 + 3z^2 + 4z + 2
    # The coefficients are ordered from the highest power to the constant term.
    coeffs = [1, 3, 4, 2]
    
    # Calculate the roots of the polynomial
    roots = np.roots(coeffs)
    
    # The polynomial p(z) is z^3 + 3z^2 + 4z + 2
    print("The Ehrhart polynomial for d=3 (scaled by a factor of 2) is:")
    print("p(z) = 1*z^3 + 3*z^2 + 4*z + 2")
    print("\nIts roots are:")
    for root in roots:
        # np.round is used to handle floating point inaccuracies
        print(np.round(root, 5))
        
    print("\nAnalyzing the properties of the roots:")
    
    # Check property A: Every root of p has real part -1.
    all_real_part_neg_one = all(np.isclose(r.real, -1) for r in roots)
    print(f"A. Do all roots have real part -1? {all_real_part_neg_one}")

    # Check property B: Every root of p is real.
    all_real = all(np.isclose(r.imag, 0) for r in roots)
    print(f"B. Are all roots real? {all_real}")
    
    # Check property E: Every root of p has real part -1/2.
    all_real_part_neg_half = all(np.isclose(r.real, -0.5) for r in roots)
    print(f"E. Do all roots have real part -1/2? {all_real_part_neg_half}")

solve_poly()