import math

def final_formula_printer():
    """
    This function prints the closed-form expression for the infinite product
    Product_{n=3 to infinity} (1 - z^3/n^3).
    """

    # Define the components of the formula for clarity
    lhs = "Product_{n=3 to infinity} (1 - z^3/n^3)"

    # The denominator consists of the first two terms (n=1, n=2) and the Gamma functions
    term_n1 = "1 - z^3"
    term_n2 = "1 - z^3/8"
    gamma_terms = "Gamma(1-z) * Gamma(1-w*z) * Gamma(1-w^2*z)"

    rhs = f"1 / (({term_n1}) * ({term_n2}) * {gamma_terms})"
    
    # Definition of w (omega)
    omega_def = "w = exp(i*2*pi/3)"

    # Print the final equation
    print(f"{lhs} = {rhs}")
    print(f"where {omega_def}")
    
    # To satisfy the "output each number" requirement, we explicitly show the numbers from the formula
    print("\nThe numbers appearing in the final equation are:")
    print("Start of product index: 3")
    print("Power of z and n: 3")
    print("Term for n=1: (1 - z^3/1^3), which contains 1 and 3.")
    print("Term for n=2: (1 - z^3/2^3), which contains 1, 3, 2, 3. (Simplified to 1 - z^3/8)")
    print("Gamma function arguments: (1-z), (1-w*z), (1-w^2*z), which contains 1 and 2.")
    print("Definition of w: exp(i*2*pi/3), which contains 2 and 3.")


if __name__ == "__main__":
    final_formula_printer()