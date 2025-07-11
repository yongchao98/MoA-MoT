def generate_equation():
    """
    Generates and prints the equation for the function f(x) based on graph analysis.
    The graph has x-intercepts at -b, d, b and vertical asymptotes at a, c.
    """
    
    # Define the symbolic variables for clarity, though they are not used in calculation
    # a, b, c, d are positive constants based on the graph
    
    # Numerator is derived from the roots: (x+b)(x-b)(x-d)
    # N(x) = (x^2 - b^2)(x - d) = x^3 - d*x^2 - b^2*x + b^2*d
    
    # Denominator is derived from the vertical asymptotes: (x-a)(x-c)
    # D(x) = x^2 - (a+c)x + a*c
    
    # The problem asks to output each number in the final equation.
    # We will print the coefficients and terms of the expanded polynomials.
    
    numerator_str = "1*x^3 - d*x^2 - b^2*x + b^2*d"
    denominator_str = "1*x^2 - (a+c)*x + a*c"
    
    # To be more explicit about the coefficients as requested:
    # Numerator coefficients: 1, -d, -b^2, b^2*d
    # Denominator coefficients: 1, -(a+c), a*c
    
    # We will print the final equation in a readable format.
    # The instruction "output each number" is interpreted as showing the full expanded equation.
    
    print(f"f(x) = ({numerator_str}) / ({denominator_str})")

generate_equation()