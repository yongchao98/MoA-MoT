def generate_equation():
    """
    This function generates and prints the equation for f(x) based on graphical analysis.
    The symbolic constants a, b, c, and d are used as placeholders for the actual numerical values.
    """
    
    # Define the symbolic constants and exponents as strings for printing
    a = "a"
    b = "b"
    c = "c"
    d = "d"
    exponent = "2"
    
    # Construct the numerator and denominator strings
    numerator = f"(x - {d})^{exponent} * (x^{exponent} - {b}^{exponent})"
    denominator = f"(x - {a})^{exponent} * (x - {c})"
    
    # Print the final equation for f(x)
    print(f"f(x) = ({numerator}) / ({denominator})")

generate_equation()