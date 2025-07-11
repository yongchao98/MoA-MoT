import numpy as np
import math

def solve_poly_properties(d):
    """
    Analyzes the Ehrhart polynomial for a given dimension d.
    """
    if d < 2:
        print("Dimension d must be 2 or greater.")
        return

    # From the derivation, the roots of the Ehrhart polynomial p(z) are
    # -1 (with multiplicity 2), -2, -3, ..., -(d-1).
    roots = [-1.0, -1.0] + [float(-i) for i in range(2, d)]
    
    # The leading coefficient is 1/(d-1)!
    leading_coeff = 1.0 / math.factorial(d - 1)
    
    # Construct the polynomial from its roots
    # np.poly returns a monic polynomial, so we scale it.
    poly_coeffs = leading_coeff * np.poly(roots)
    
    # Create a polynomial object for easier handling
    p = np.poly1d(poly_coeffs)
    
    print(f"For dimension d = {d}:")
    print("-" * 25)

    # Print the polynomial equation
    # We will format the polynomial string manually for clarity
    poly_string = "p(z) = "
    for i, coeff in enumerate(p.coeffs):
        power = p.order - i
        if abs(coeff) > 1e-9: # a tolerance for floating point
            sign = " - " if coeff < 0 else " + "
            if i == 0:
                sign = "" if coeff > 0 else "-"
            
            val = abs(coeff)
            
            if power > 1:
                poly_string += f"{sign}{val:.4f}*z^{power}"
            elif power == 1:
                poly_string += f"{sign}{val:.4f}*z"
            else:
                poly_string += f"{sign}{val:.4f}"
    print(poly_string)
    
    # The final equation requested by the prompt can be the factored form:
    print("\nThe factored form of the polynomial is:")
    factors = ""
    for r in sorted(list(set(roots))):
        multiplicity = roots.count(r)
        if multiplicity > 1:
            factors += f"(z + {abs(r)})^{multiplicity}"
        else:
            factors += f"(z + {abs(r)})"
    print(f"p(z) = (1/{math.factorial(d-1)}) * {factors}")


    print("\nRoots of the polynomial:", np.roots(p))
    print("-" * 25)
    
    # Check the given options
    print("Checking the answer choices:")
    
    # A. Every root of p has real part -1.
    is_A = all(np.real(r) == -1 for r in roots)
    print(f"A. Every root of p has real part -1.  --> {is_A}")

    # B. Every root of p is real.
    is_B = all(np.isreal(r) for r in roots)
    print(f"B. Every root of p is real.          --> {is_B}")

    # C. The coefficients of p sum exactly d.
    # The sum of coefficients is p(1)
    sum_coeffs = p(1)
    is_C = math.isclose(sum_coeffs, d)
    print(f"C. The sum of coefficients is {sum_coeffs:.1f}, not {d}. --> {is_C}")

    # D. The coefficients of p sum exactly d choose d/2.
    if d % 2 == 0:
        d_choose_d_half = math.comb(d, d // 2)
        is_D = math.isclose(sum_coeffs, d_choose_d_half)
        print(f"D. The sum of coefficients is {sum_coeffs:.1f}, not {d_choose_d_half}. --> {is_D}")
    else:
        is_D = False
        print(f"D. d is odd, so d/2 is not an integer. --> {is_D}")

    # E. Every root of p has real part -1/2.
    is_E = all(np.real(r) == -0.5 for r in roots)
    print(f"E. Every root of p has real part -1/2. --> {is_E}")

# Let's run the analysis for d=4 as an example
solve_poly_properties(d=4)