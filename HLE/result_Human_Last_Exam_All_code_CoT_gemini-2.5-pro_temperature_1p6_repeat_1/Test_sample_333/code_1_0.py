import sys

def solve_projectile_problem():
    """
    Calculates the landing position of the second fragment of an exploded projectile.
    """
    # Given horizontal distance to the highest point of elevation in meters.
    I = 500

    # The problem can be solved using the principle of the conservation of momentum.
    # The center of mass (CM) of the system of two fragments continues on the
    # original trajectory of the projectile, as the explosion is an internal force.
    #
    # 1. The total range (R) of the unexploded projectile would be 2 * I.
    #    This is where the CM of the two fragments will land.
    #    R = 2 * I
    #
    # 2. The projectile explodes into two equal parts. The landing position of
    #    their CM is the average of their individual landing positions (x1 and x2).
    #    R = (x1 + x2) / 2
    #
    # 3. We are told one fragment fell near the gun, so its landing position x1 = 0.
    #
    # 4. Substituting the values into the equation:
    #    2 * I = (0 + x2) / 2
    #    This gives us the final formula: x2 = 4 * I

    # Calculate the landing position of the second fragment.
    multiplier = 4
    x2 = multiplier * I

    # Print the explanation and the final calculation.
    print(f"The landing position of the second fragment (x_2) can be determined from the motion of the center of mass.")
    print(f"The formula is: x_2 = 4 * I")
    print(f"Given the distance to the highest point, I = {I} m.")
    print(f"The final calculation is:")
    # Using sys.stdout.write to prevent the print function from adding a space after the numbers
    sys.stdout.write(f"{multiplier} * {I} = {x2}\n")
    print(f"\nTherefore, the maximum distance from the gun to be safe from the second fragment is {x2} meters.")

solve_projectile_problem()

# The final answer is the numerical result of the calculation.
print("<<<2000>>>")