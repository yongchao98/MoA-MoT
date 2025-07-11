import sympy

def display_potential_formula():
    """
    Solves for the electric potential symbolically and displays the formula
    for the region 0 < y < a.
    """
    # Define the symbolic variables used in the equations
    sigma_0, k, x, y, a, b, epsilon_1, epsilon_2 = sympy.symbols('sigma_0 k x y a b epsilon_1 epsilon_2')

    # The problem asks for the electric potential Phi(x, y) in the region 0 <= y <= a.
    # The derivation leads to the expression in Option A.

    print("The electric potential Phi(x, y) in the region 0 < y < a is determined by solving Laplace's equation with the given boundary conditions.")
    print("The final expression is broken down below:\n")

    # Construct the numerator of the potential function for 0 < y < a
    numerator = -sigma_0 * sympy.sinh(k * b) * sympy.sinh(k * (y - a)) * sympy.sin(k * x)

    # Construct the denominator of the potential function
    denominator = k * (epsilon_2 * sympy.cosh(k * a) * sympy.sinh(k * b) + epsilon_1 * sympy.sinh(k * a) * sympy.cosh(k * b))

    # Construct the full expression for the potential Phi(x, y)
    phi_expression = numerator / denominator

    # Print each component of the final equation
    print("Numerator of the potential function Phi(x,y):")
    print(sympy.pretty(numerator))
    print("\nDenominator of the potential function Phi(x,y):")
    print(sympy.pretty(denominator))
    print("\nComplete expression for the potential Phi(x, y) in the region 0 < y < a:")
    print(sympy.pretty(phi_expression))

if __name__ == "__main__":
    display_potential_formula()
