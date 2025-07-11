import sympy

def solve_minimal_polynomial():
    """
    This function calculates the minimal polynomial for the shortest geodesic
    distance an ant can walk on a dodecahedron from a vertex back to itself.

    The steps are:
    1. Define the shortest distance L = (sqrt(5) - 1) / 2 based on geometric derivation.
    2. Use the sympy library to compute the minimal polynomial for L.
    3. Print the resulting polynomial equation, showing each coefficient as requested.
    """
    
    # The variable of the polynomial is 'x'
    x = sympy.Symbol('x')

    # The shortest distance L, as derived from the geometry of the unfolded dodecahedron,
    # is (sqrt(5) - 1) / 2.
    shortest_distance = (sympy.sqrt(5) - 1) / 2

    # sympy's minimal_polynomial function finds the polynomial of lowest degree
    # with rational (in this case, integer) coefficients that has the given number as a root.
    min_poly = sympy.minimal_polynomial(shortest_distance, x)

    # Extract the coefficients and the degree of the polynomial for printing.
    coeffs = min_poly.all_coeffs()
    degree = min_poly.degree()
    
    print("The minimal polynomial P(x) = 0 for the shortest distance x is:")
    
    output_parts = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        # Skip any term with a zero coefficient
        if coeff == 0:
            continue
        
        # Determine the sign of the term
        # The first term's sign is part of the coefficient
        if i > 0:
            if coeff > 0:
                output_parts.append("+")
            else:
                output_parts.append("-")
        elif coeff < 0:
             output_parts.append("-")
            
        # Add the absolute value of the coefficient
        # The problem asks for each number in the equation.
        output_parts.append(str(abs(int(coeff))))
        
        # Add the variable 'x' and its exponent if it's not a constant term
        if power > 0:
            output_parts.append("*x")
            if power > 1:
                output_parts.append(f"**{power}")

    # Join all parts with spaces and print the final equation
    print(" ".join(output_parts), "= 0")

solve_minimal_polynomial()