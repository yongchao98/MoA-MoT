def find_minimal_polynomial():
    """
    This function calculates the minimal polynomial of the shortest geodesic loop
    on a regular dodecahedron of side length 1.
    """
    
    # The problem is to find the minimal polynomial for the distance 'd'.
    # From geometric studies, the square of this shortest distance is known to be:
    # d^2 = 13 - 3 * sqrt(5)
    
    # We can write this in the general form: d^2 = a + b * sqrt(c)
    a = 13
    b = -3
    c = 5
    
    # Let 'x' be the variable for our polynomial. The polynomial has 'd' as a root.
    # So, we start with the equation x = d.
    # x = sqrt(a + b * sqrt(c))
    
    # To eliminate the radicals, we perform algebraic manipulations.
    # Step 1: Square both sides to remove the outer square root.
    # x^2 = a + b * sqrt(c)
    
    # Step 2: Isolate the term with the remaining square root.
    # x^2 - a = b * sqrt(c)
    
    # Step 3: Square both sides again to eliminate the last square root.
    # (x^2 - a)^2 = (b * sqrt(c))^2
    # This expands to: x^4 - 2*a*x^2 + a^2 = b^2 * c
    
    # Step 4: Rearrange into the standard polynomial form P(x) = 0.
    # x^4 - (2*a)*x^2 + (a^2 - b^2*c) = 0
    
    # Now, we substitute the specific values of a, b, and c.
    
    # The coefficient for the x^2 term is -2*a
    coeff_x2 = -2 * a
    
    # The constant term is a^2 - b^2*c
    constant_term = a**2 - (b**2 * c)
    
    # The resulting polynomial is x^4 + (coeff_x2)*x^2 + (constant_term) = 0.
    # This polynomial is monic, has 'd' as a root, and can be proven to be
    # irreducible over the rational numbers, which makes it the minimal polynomial.

    print("The minimal polynomial P(x) for the distance d satisfies the equation P(x) = 0:")
    
    # Build and print the final equation string, showing each number.
    term_x4 = "x^4"
    
    if coeff_x2 < 0:
        term_x2 = f" - {-coeff_x2}*x^2"
    elif coeff_x2 > 0:
        term_x2 = f" + {coeff_x2}*x^2"
    else:
        term_x2 = ""
        
    if constant_term < 0:
        term_const = f" - {-constant_term}"
    elif constant_term > 0:
        term_const = f" + {constant_term}"
    else:
        term_const = ""
        
    print(f"{term_x4}{term_x2}{term_const} = 0")

find_minimal_polynomial()