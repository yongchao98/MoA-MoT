import sys

def solve_projectile_explosion():
    """
    Calculates the landing position of the second fragment of an exploded projectile.
    """
    # I is the horizontal distance from the gun to the highest point of elevation in meters.
    I = 500

    # One fragment fell near the gun. We take its landing position to be 0.
    x1 = 0

    # The explosion is from internal forces, so the Center of Mass (CM) continues
    # on the original trajectory. The projectile's full range would be 2 * I.
    # This is the landing spot of the center of mass.
    x_cm = 2 * I

    # The formula for the center of mass of two equal fragments is:
    # x_cm = (x1 + x2) / 2
    # We need to find x2, the landing spot of the second fragment.
    # Rearranging the formula: x2 = 2 * x_cm - x1
    x2 = 2 * x_cm - x1

    # Print the final equation with all numbers included
    print("The final safe distance (x2) is calculated using the center of mass principle.")
    print("The landing position of the Center of Mass (x_cm) is twice the distance to the peak:")
    print(f"x_cm = 2 * {I} = {x_cm} m")
    print("\nThe landing position of the second fragment (x2) is then found with the formula:")
    print("x2 = 2 * x_cm - x1")
    print("Substituting the values:")
    print(f"x2 = 2 * {x_cm} - {x1}")
    print(f"x2 = {x2} m")
    print("\nTherefore, the maximum distance from the gun you can safely be is 2000 meters.")
    
    # Writing the final answer in the requested format to stdout
    # This part is for automated checking and might not be visible in all environments.
    sys.stdout.write(f"\n<<<{int(x2)}>>>")

solve_projectile_explosion()