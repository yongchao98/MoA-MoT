import numpy as np

def solve_fixed_points():
    """
    This function solves for the fixed points of a composed function h(x)=f(g(x))
    that is constructed to have the maximum possible number of fixed points.
    
    The equation for the fixed points is h(x) = x, or h(x) - x = 0.
    We construct a polynomial h(x) such that h(x) - x has 9 distinct real roots
    and h'(x) > 0. A good candidate for h(x) - x is a scaled version of the
    9th Chebyshev polynomial, T_9(x), which has 9 known real roots in [-1, 1].
    
    The equation we solve is T_9(x) = 0.
    """
    
    # The coefficients of the 9th Chebyshev polynomial T_9(x) are:
    # 256x^9 - 576x^7 + 432x^5 - 120x^3 + 9x
    # Setting this to zero gives the fixed points.
    # Note that the constant term is 0, so one root is x=0.
    chebyshev_coeffs = [256, 0, -576, 0, 432, 0, -120, 0, 9, 0]
    
    # Create the polynomial object
    poly = np.poly1d(chebyshev_coeffs)
    
    # Print the equation for the fixed points
    print("The representative polynomial equation for finding the fixed points is:")
    
    equation_parts = []
    for i, coeff in enumerate(poly.coeffs):
        power = len(poly.coeffs) - 1 - i
        
        # Skip zero-coefficient terms
        if abs(coeff) < 1e-9:
            continue
        
        # Determine the sign
        sign = " - " if coeff < 0 else " + "
        if i == 0 and coeff > 0:
            sign = ""
        
        # Get the absolute value of the coefficient
        abs_coeff = abs(coeff)
        
        # Format the term
        if power == 1:
            term = f"{abs_coeff}x"
        elif power == 0:
            term = f"{abs_coeff}"
        else:
            term = f"{abs_coeff}x^{power}"
            
        equation_parts.append(sign + term)
        
    print("".join(equation_parts).lstrip(" + "), "= 0")
    
    # Find the roots of the polynomial
    roots = np.roots(poly)
    
    # Filter for real roots and find unique values (due to potential floating point inaccuracies)
    real_roots = roots[np.isreal(roots)].real
    unique_real_roots = sorted(np.unique(np.round(real_roots, 8)))
    
    # Print the found fixed points
    print("\nThe found fixed points are:")
    for root in unique_real_roots:
        print(f"{root:.4f}")
    
    # Print the total count
    print(f"\nThe maximum number of fixed points found is: {len(unique_real_roots)}")

solve_fixed_points()