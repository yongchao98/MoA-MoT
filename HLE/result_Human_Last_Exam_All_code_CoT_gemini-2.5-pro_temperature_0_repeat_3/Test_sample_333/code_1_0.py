def solve_projectile_explosion():
    """
    Calculates the maximum landing distance of a fragment from an exploding projectile.
    """
    # Given horizontal distance to the highest point of elevation
    I = 500  # in meters

    # The total range (R) the projectile would have traveled without exploding is twice the distance to the highest point.
    # This is also the landing position of the center of mass of the fragments.
    R = 2 * I

    # One fragment fell near the gun. To maximize the distance of the second fragment,
    # we assume the first fragment lands at the gun's position.
    x1 = 0  # in meters

    # The position of the center of mass (R) for two equal-mass fragments is the average of their landing positions (x1 and x2).
    # The equation is: R = (x1 + x2) / 2
    # We solve for x2: x2 = 2 * R - x1
    x2 = 2 * R - x1

    print(f"The horizontal distance to the highest point is I = {I} m.")
    print(f"The total range the center of mass will travel is R = 2 * I = {R} m.")
    print(f"The first fragment is assumed to land at x1 = {x1} m to maximize the second fragment's range.")
    print("\nThe governing equation is: R = (x1 + x2) / 2")
    print("Substituting the values to show the final calculation:")
    # We use int() to display whole numbers for clarity in the final equation
    print(f"{R} = ({int(x1)} + {int(x2)}) / 2")
    print(f"\nTherefore, the maximum distance from the gun the second fragment can land is {int(x2)} m.")

solve_projectile_explosion()