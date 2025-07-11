import math

def generate_formula_string():
    """
    This function generates and prints the symbolic formula for the infinite product.
    The derivation involves the Gamma function and adjusting the product's starting index.
    """
    
    # The left-hand side of the equation is the product itself.
    # The numbers 3 and 3 are explicitly written.
    lhs = "Product_{n=3 to infinity} (1 - z^3/n^3)"
    
    # The term for n=1 is (1 - z^3/1^3). The number 1 and 3 are present.
    # The term for n=2 is (1 - z^3/2^3). 2^3 = 8. The numbers 2 and 3 -> 8 are here.
    # The simplification leads to 8 in the numerator.
    rhs_numerator = "8"
    
    # The denominator consists of the polynomial part from the n=1 and n=2 terms
    # and the Gamma function part from the product starting at n=1.
    denominator = "(1 - z^3) * (8 - z^3) * Gamma(1 - z) * Gamma(1 - z*w) * Gamma(1 - z*w^2)"
    
    # Definition of w (omega) which contains numbers 2 and 3.
    w_definition = "where w = exp(2*pi*i/3) is a primitive cube root of unity."
    
    # Construct the final string representation of the equation.
    equation = f"{lhs} = {rhs_numerator} / ({denominator})"
    
    full_answer = f"{equation},\n{w_definition}"
    
    print(full_answer)

generate_formula_string()
