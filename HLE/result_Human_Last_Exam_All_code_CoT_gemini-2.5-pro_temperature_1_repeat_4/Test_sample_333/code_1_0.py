def solve_projectile_explosion():
    """
    Calculates the maximum distance a second fragment can travel after a projectile
    explodes at its highest point.
    """
    # Given horizontal distance to the highest point of elevation in meters.
    I = 500

    # According to the principle of the conservation of the center of mass,
    # the landing position of the second fragment (D) can be determined by the formula:
    # D = 4 * I
    # where I is the distance to the peak and one fragment lands at the origin (0).

    # The constant multiplier in our derived formula.
    multiplier = 4

    # Calculate the maximum distance for the second fragment.
    max_distance = multiplier * I

    # Print the explanation and the final equation with all its components.
    print("The problem is solved using the center of mass principle.")
    print("The formula for the maximum distance (D) is D = 4 * I.")
    print("The final equation with the given values is:")
    print(f"{max_distance} m = {multiplier} * {I} m")
    print("\nTherefore, the maximum safe distance from the gun is:")
    print(f"{max_distance} meters.")

solve_projectile_explosion()