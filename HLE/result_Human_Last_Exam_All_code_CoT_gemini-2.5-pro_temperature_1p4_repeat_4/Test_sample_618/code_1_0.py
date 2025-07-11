import sympy

def display_infinitesimal_transformation():
    """
    This function calculates and displays the general representation for the
    infinitesimal transformation on the variable 'x' for the given heat equation.
    """
    # Define symbols for the constants and variables
    # c_2 and c_4 are arbitrary constants derived from the symmetry analysis.
    # k_1 is the constant parameter from the PDE.
    # t is the time variable.
    c_2, c_4, k_1, t = sympy.symbols('c_2 c_4 k_1 t')

    # The general representation for the infinitesimal transformation on x, denoted as xi.
    # This result is obtained by solving the determining equations from Lie symmetry analysis.
    xi = c_4 - (2 * c_2 / k_1) * sympy.exp(k_1 * t)

    print("The general representation for the infinitesimal transformation on x, xi(t), is:")
    # Use sympy.pretty_print for a clear mathematical representation
    sympy.pretty_print(xi)

    print("\nWhere:")
    print(f"  - c_2 and c_4 are arbitrary constants representing different symmetries.")
    print(f"  - k_1 is the non-zero constant from the differential equation's source term.")
    print(f"  - t is the time variable.")

    # Fulfilling the requirement to output each number in the final equation.
    # We parse the expression to find numerical constants.
    # The only explicit number in the derived formula is '2'.
    # Note: the coefficient of c_4 is implicitly 1, and the exponent of exp is 1.
    numbers = [s for s in xi.atoms(sympy.Number)]
    print("\nThe explicit numerical constant in the final equation is:")
    for num in numbers:
        # We only care about the '2' which is explicitly part of the structure
        if num == 2:
            print(int(num))

if __name__ == '__main__':
    display_infinitesimal_transformation()
