import math

def solve_projectile_problem():
    """
    Calculates the landing position of the second fragment of an exploded projectile.
    """
    # I: Horizontal distance from the gun to the highest point of elevation (m)
    I = 500

    # x1: Landing position of the first fragment (m)
    # It fell "near the gun", so we assume its landing position is at the origin.
    x1 = 0

    # The center of mass of the system continues on the original trajectory.
    # The full range (R) of the original trajectory is twice the distance to the highest point.
    # This is also the final landing position of the center of mass.
    R_cm = 2 * I

    # For two equal masses, the center of mass is the average of their positions:
    # R_cm = (x1 + x2) / 2
    # We solve for x2, the landing position of the second fragment.
    # x2 = 2 * R_cm - x1
    x2 = 2 * R_cm - x1

    print("This script calculates the maximum safe distance from a gun after a projectile explosion.")
    print("The calculation is based on the principle of conservation of the center of mass trajectory.")
    print(f"\nGiven values:")
    print(f"Horizontal distance to highest point (I): {I} m")
    print(f"Landing position of the first fragment (x1): {x1} m")
    
    print("\nStep 1: Calculate the range of the center of mass (R_cm).")
    print(f"R_cm = 2 * I = 2 * {I} = {R_cm}")

    print("\nStep 2: Use the center of mass formula R_cm = (x1 + x2) / 2 to find x2.")
    print("The final equation is: x2 = 2 * R_cm - x1")
    print(f"Final Calculation: {2} * {R_cm} - {x1} = {int(x2)}")
    
    print(f"\nTherefore, the second fragment will land at {int(x2)} m from the gun.")
    print(f"The maximum safe distance is {int(x2)} m.")

solve_projectile_problem()