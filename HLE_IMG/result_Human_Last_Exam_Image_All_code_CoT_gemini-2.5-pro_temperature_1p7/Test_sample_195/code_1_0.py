def generate_equation():
    """
    This function generates and prints the equation for f(x) based on the graphical analysis.
    The parameters a, b, and c from the graph define the function.
    The root 'd' is determined by the relation d = a + c.
    """

    # Define the parameters as symbolic strings for the output
    a, b, c = 'a', 'b', 'c'

    # The equation is expressed in a factored form to clearly show the roots and poles.
    # The requirement to "output each number in the final equation" is interpreted as
    # showing each parameter that defines the function's characteristic points (roots and poles).

    numerator = f"(x - {b})(x + {b})(x - ({a} + {c}))"
    denominator = f"(x - {a})(x - {c})"
    
    equation = f"f(x) = {numerator} / {denominator}"

    print(equation)

generate_equation()