import math

def print_infinite_product_formula():
    """
    This function prints the closed-form expression for the infinite product:
    Product_{n=3 to infinity} (1 - z^3 / n^3)

    The formula is derived using the infinite product representation of the Gamma function.
    """
    
    # Define the components of the formula
    numerator = "1"
    
    # Term for n=1
    term_n1 = "(1 - z**3 / 1**3)"
    
    # Term for n=2
    term_n2 = "(1 - z**3 / 2**3)"
    
    # Gamma function part from the product starting at n=1
    # where omega = exp(2*pi*i/3) and omega^2 = exp(4*pi*i/3)
    gamma_product = "Gamma(1 - z) * Gamma(1 - z*exp(2*pi*i/3)) * Gamma(1 - z*exp(4*pi*i/3))"
    
    # The full denominator is the product of these terms
    denominator = f"{term_n1} * {term_n2} * {gamma_product}"

    # Print the final result in a structured way
    print("The closed-form expression for the infinite product is:")
    print()
    print(f"  {numerator}")
    print("---------------------------------------------------------------------------------------------------")
    print(f"  ({term_n1}) * ({term_n2}) * {gamma_product}")
    print()
    print("Which can be written as:")
    print("1 / ( " + denominator + " )")
    print()
    print("Where:")
    print("  - z is a complex variable.")
    print("  - Gamma is the Gamma function.")
    print("  - exp is the exponential function (e^x).")
    print("  - pi is the mathematical constant Pi (e.g., math.pi).")
    print("  - i is the imaginary unit (the square root of -1).")

if __name__ == '__main__':
    print_infinite_product_formula()