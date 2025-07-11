def solve_projectile_explosion():
    """
    Calculates the landing position of the second fragment of an exploded projectile.

    The solution is based on the principle that the center of mass of the system
    continues on its original trajectory unaffected by the internal forces of the explosion.
    """
    # Horizontal distance from the gun to the highest point of elevation in meters.
    I = 500

    # The landing position of the first fragment, which falls near the gun.
    x1 = 0

    # The center of mass (CM) of the system continues on the original path.
    # The total range of the original projectile would be 2 * I.
    # So, the landing position of the center of mass is x_cm = 2 * I.
    x_cm = 2 * I

    # The formula for the center of mass for two equal fragments is x_cm = (x1 + x2) / 2.
    # We can solve for x2, the landing position of the second fragment:
    # 2 * x_cm = x1 + x2
    # x2 = 2 * x_cm - x1
    # Substituting x_cm = 2 * I:
    # x2 = 2 * (2 * I) - x1 = 4 * I - x1

    # Calculate the final distance for the second fragment.
    x2 = 4 * I - x1

    print("The problem is solved using the concept of the center of mass.")
    print("The landing position of the second fragment (x2) is determined by the equation: x2 = 4 * I - x1")
    print(f"Where I (distance to highest point) = {I} m and x1 (landing spot of first fragment) = {x1} m.")
    print("\nFinal Equation:")
    print(f"x2 = 4 * {I} - {x1}")
    print(f"x2 = {4 * I}")

    print(f"\nThe maximum safe distance from the gun is {x2} m.")

solve_projectile_explosion()
<<<2000>>>