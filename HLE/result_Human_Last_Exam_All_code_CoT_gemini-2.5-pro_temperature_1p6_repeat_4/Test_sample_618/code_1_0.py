import sympy

def find_x_transformation():
    """
    This function uses the sympy library to define and print the general form of the 
    infinitesimal transformation on the x-coordinate for the given PDE.
    """
    
    # Define the symbolic variables and constants
    # t: time variable
    # k1: constant parameter from the PDE
    # C1, C2: arbitrary constants from solving the determining equations
    t, k1 = sympy.symbols('t k1')
    C1, C2 = sympy.symbols('C1 C2', real=True)
    
    # The general representation for the infinitesimal xi(t), which corresponds to the
    # transformation on the x-coordinate.
    # It consists of two parts:
    # 1. C1: Represents a constant spatial translation (symmetry of the PDE).
    # 2. C2 * exp(k1 * t): Represents a time-dependent, exponential scaling/translation.
    xi = C1 + C2 * sympy.exp(k1 * t)

    # Creating a symbolic representation of the function xi(t)
    xi_function = sympy.Function('ξ')(t)

    # Formatting the final equation to be printed.
    # The left side is ξ(t), representing the infinitesimal for x.
    # The right side is the derived expression.
    final_equation = sympy.Eq(xi_function, xi)

    # Print the final result in a clear format
    print("The general representation for the infinitesimal transformation on x is given by the function ξ(t):")
    
    # The sympy.pretty_print function provides a more readable output format.
    sympy.pretty_print(final_equation)
    
find_x_transformation()