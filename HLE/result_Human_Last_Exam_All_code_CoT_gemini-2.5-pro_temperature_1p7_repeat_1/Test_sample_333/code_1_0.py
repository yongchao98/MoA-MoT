def solve_projectile_explosion():
    """
    Calculates the maximum landing distance of the second fragment
    of an exploded projectile.
    """
    # The horizontal distance from the gun to the highest point of elevation.
    I = 500  # in meters

    # --- Step 1: Determine the landing spot of the Center of Mass (CM) ---
    # The CM follows the original trajectory, landing at twice the distance to the apex.
    R_cm = 2 * I

    # --- Step 2: Determine the landing spot of the first fragment (x1) ---
    # To maximize the distance of the second fragment, we assume the first fragment
    # lands at the minimum possible distance, which is back at the gun.
    x1 = 0  # in meters

    # --- Step 3: Calculate the landing spot of the second fragment (x2) ---
    # For two equal masses, R_cm = (x1 + x2) / 2.
    # Therefore, x2 = 2 * R_cm - x1.
    # Substituting R_cm = 2 * I, we get x2 = 2 * (2 * I) - x1 = 4 * I - x1.
    x2 = 4 * I - x1

    # --- Step 4: Print the final result and the equation ---
    print("The maximum distance from the gun where the second fragment can land is calculated as follows:")
    # The user requested to see the numbers in the final equation.
    # The most direct relationship is x2 = 4 * I
    print(f"Maximum Distance = 4 * I")
    print(f"Maximum Distance = 4 * {I} = {x2} meters")

solve_projectile_explosion()