def generate_equation():
    """
    This function prints the derived equation for the graph f(x).
    The variables a, b, c, d are treated as symbolic placeholders.
    """
    numerator_str = "-x^3 + d*x^2 + b^2*x - b^2*d"
    denominator_str = "x^2 - (a + c)*x + a*c"
    
    print(f"f(x) = ({numerator_str}) / ({denominator_str})")

generate_equation()