def solve_projectile_explosion():
    """
    Calculates the landing distance of the second fragment of an exploded projectile.
    """
    # Initial horizontal distance to the highest point (apex) in meters.
    I = 500

    # The total distance from the gun to the landing point of the second fragment
    # can be expressed as a multiple of I.
    # As derived from the principles of conservation of momentum:
    # Let D_total be the final distance.
    # D_total = I (distance to explosion) + 3 * I (distance after explosion)
    # D_total = 4 * I
    multiplier = 4
    
    # Calculate the final distance.
    D_total = multiplier * I

    print("The problem is solved using the principle of conservation of momentum.")
    print("Let I be the horizontal distance from the gun to the point of explosion (the apex).")
    print(f"Given I = {I} m.")
    print("One fragment falls near the gun, which means its horizontal velocity after the explosion is reversed.")
    print("By conserving momentum, the second fragment's horizontal velocity becomes 3 times the original.")
    print("The total distance (D_total) traveled by the second fragment is the distance to the apex plus the distance it travels after the explosion.")
    print("The final calculation is based on the equation: D_total = 4 * I")
    print("\nFinal Equation:")
    print(f"{D_total} m = {multiplier} * {I} m")


solve_projectile_explosion()
<<<2000>>>