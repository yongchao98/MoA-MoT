import math

def generate_expression_string():
    """
    This function generates and prints the symbolic expression for the given infinite product.
    """
    # The product we want to evaluate is P = Prod_{n=3 to infinity} (1 - z^3 / n^3)
    # This can be expressed using the Gamma function.

    # The formula for the product from n=1 is:
    # P_1 = 1 / (Gamma(1 - z) * Gamma(1 - z*w) * Gamma(1 - z*w^2))
    # where w = exp(i*2*pi/3)

    # We need to divide P_1 by the terms for n=1 and n=2.
    # Term n=1: (1 - z^3/1^3)
    # Term n=2: (1 - z^3/2^3)

    # The final expression is therefore:
    lhs_str = "Product_{n=3 to infinity} (1 - z^3/n^3)"
    
    # Using 'w' for exp(i*2*pi/3) and 'w^2' for exp(i*4*pi/3) for clarity
    # rhs_str = "1 / ((1 - z^3) * (1 - z^3/8) * Gamma(1 - z) * Gamma(1 - z*w) * Gamma(1 - z*w^2))"
    # To be more explicit with numbers as requested:
    rhs_str = "1 / ((1 - z^3) * (1 - z^3/8) * Gamma(1 - z) * Gamma(1 - z*exp(i*2*pi/3)) * Gamma(1 - z*exp(i*4*pi/3)))"

    # We output the full equation to show the result clearly.
    # All numbers from the equation (3, 3, 1, 8, 1, 2, 3, 4) are present in the output string.
    final_equation = f"{lhs_str} = {rhs_str}"
    print(final_equation)

generate_expression_string()