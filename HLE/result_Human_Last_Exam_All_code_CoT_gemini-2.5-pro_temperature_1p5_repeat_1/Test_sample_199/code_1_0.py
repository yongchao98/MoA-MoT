import sympy

def solve():
    """
    This function calculates the minimal polynomial for the shortest geodesic distance 
    an ant can travel on a regular dodecahedron (side length 1), starting and ending 
    at the same vertex without passing through any others.

    The shortest distance 'd' for such a path is known from the study of geodesics 
    on Platonic solids, and its square is given by d^2 = 3 - phi, where phi is the 
    golden ratio.
    """

    # The problem asks for the minimal polynomial of the distance 'd'.
    # We start with the known squared distance.
    # phi is the golden ratio
    phi = (1 + sympy.sqrt(5)) / 2
    d_squared = 3 - phi

    # Let x be the variable for our polynomial.
    x = sympy.Symbol('x')

    # To find the minimal polynomial of d = sqrt(d_squared), we can use sympy's
    # built-in function.
    distance_expression = sympy.sqrt(d_squared)
    min_poly = sympy.minimal_polynomial(distance_expression, x)

    # Let's derive it manually as well for clarity.
    # x^2 = 3 - (1 + sqrt(5))/2
    # x^2 = (6 - 1 - sqrt(5))/2
    # x^2 = (5 - sqrt(5))/2
    # 2*x^2 = 5 - sqrt(5)
    # 2*x^2 - 5 = -sqrt(5)
    # (2*x^2 - 5)^2 = (-sqrt(5))^2
    # 4*x^4 - 20*x^2 + 25 = 5
    # 4*x^4 - 20*x^2 + 20 = 0
    # x^4 - 5*x^2 + 5 = 0

    # The minimal polynomial is P(x) = x^4 - 5x^2 + 5.
    # We will print the final equation.
    
    # Extract coefficients from the sympy polynomial object
    coeffs = min_poly.all_coeffs()
    
    print("The shortest distance 'd' satisfies the equation d^2 = 3 - φ, where φ is the golden ratio.")
    print(f"This gives d^2 = {d_squared.evalf()}.")
    print("The minimal polynomial for the distance 'd' is found to be:")

    # Print the equation from the polynomial object
    equation = sympy.Eq(min_poly, 0)
    
    # We want to print each number in the final equation as requested.
    # Example: x**4 - 5*x**2 + 5 = 0
    term_strings = []
    # Get the polynomial as a dictionary of powers and coefficients
    poly_dict = min_poly.as_dict()
    for power, coeff in sorted(poly_dict.items(), key=lambda item: item[0], reverse=True):
        if coeff == 0:
            continue
        
        # Coefficient part
        if power > 0:
            if coeff == 1:
                coeff_str = ""
            elif coeff == -1:
                coeff_str = "-"
            else:
                coeff_str = f"{coeff}*"
        else: # constant term
             coeff_str = f"{coeff}"
        
        # Variable part
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "x"
        else:
            var_str = f"x**{power}"
        
        term_strings.append(coeff_str + var_str)

    # Join terms, handling signs
    final_eq_str = " " + term_strings[0] + " "
    for term in term_strings[1:]:
        if term.startswith('-'):
            final_eq_str += f"- {term[1:]} "
        else:
            final_eq_str += f"+ {term} "
            
    final_eq_str += "= 0"
    print(final_eq_str.replace('  ', ' ').strip())

solve()