import sympy

def solve_minimal_polynomial():
    """
    This function calculates and prints the minimal polynomial for the shortest
    geodesic distance on a dodecahedron as described in the problem.
    """
    
    # Step 1: Define the problem's variable.
    # The shortest non-trivial geodesic loop from a vertex back to itself on a
    # regular dodecahedron with side length 1 has a distance 'd'. The square
    # of this distance is known from geometry to be d^2 = 7 - 2*sqrt(5).
    # We need to find the minimal polynomial of d.

    # Step 2: Use SymPy for symbolic manipulation.
    # Let 'x' be the variable of our polynomial.
    x = sympy.Symbol('x')

    # Define the value of the distance 'd' symbolically.
    d_val = sympy.sqrt(7 - 2 * sympy.sqrt(5))

    # Step 3: Let the distance 'd' be a root 'x' of the polynomial and derive it.
    # x = sqrt(7 - 2*sqrt(5))
    # x^2 = 7 - 2*sqrt(5)
    # x^2 - 7 = -2*sqrt(5)
    # Squaring both sides to eliminate the square root:
    # (x^2 - 7)^2 = (-2*sqrt(5))^2
    # x^4 - 14*x^2 + 49 = 4 * 5
    # x^4 - 14*x^2 + 49 = 20
    # x^4 - 14*x^2 + 29 = 0
    # This is the minimal polynomial. We can verify this with SymPy's built-in function.
    
    min_poly = sympy.minimal_polynomial(d_val, x)

    # Step 4: Format and print the final equation.
    # The prompt requires that each number in the final equation be outputted.
    # We will build the equation string from the polynomial's coefficients.
    
    coeffs = min_poly.all_coeffs() # Returns coefficients from highest degree to lowest
    degree = min_poly.degree()
    
    equation_parts = []
    
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        # Skip terms with a zero coefficient
        if coeff == 0:
            continue
            
        # Determine the sign
        if i > 0 and coeff > 0:
            sign = " + "
        elif coeff < 0:
            sign = " - "
        else: # First term
            sign = ""
            
        num = abs(coeff)
        
        # Format the variable part (e.g., x^4, x)
        if power > 1:
            var_part = f"x^{power}"
        elif power == 1:
            var_part = "x"
        else: # power is 0
            var_part = ""

        # Construct the term string, ensuring the coefficient number is always present
        if var_part:
            term = f"{sign}{num}*{var_part}"
        else:
            term = f"{sign}{num}"

        equation_parts.append(term)
        
    final_equation = "".join(equation_parts) + " = 0"
    
    print("The minimal polynomial for the shortest possible distance is:")
    print(final_equation)

if __name__ == "__main__":
    solve_minimal_polynomial()
