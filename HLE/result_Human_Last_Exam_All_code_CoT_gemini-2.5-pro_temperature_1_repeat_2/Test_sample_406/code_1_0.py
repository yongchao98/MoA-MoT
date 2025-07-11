def print_tripled_fixed_point_conditions():
    """
    Prints the defining equations for an FGH-tripled fixed point.
    
    A point (x, y, z) is an FGH-tripled fixed point if it satisfies
    the following system of equations, based on the standard cyclic
    definition of the functions F, G, and H.
    """
    
    # The list of conditions (equations) that must be met.
    conditions = [
        "F(x, y, z) = x",
        "G(y, z, x) = y",
        "H(z, x, y) = z"
    ]
    
    print("The conditions for a point (x, y, z) to be an FGH-tripled fixed point are given by the following system of equations:")
    
    # We iterate through the list and print each equation, prefixed with its number.
    for i, eq in enumerate(conditions):
        equation_number = i + 1
        print(f"Equation {equation_number}: {eq}")

# Execute the function to display the conditions.
print_tripled_fixed_point_conditions()