import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # Step 1: Define the parameters of the two metric spaces.
    # Space X: The interval [0, 1]. Its length is L.
    L = 1.0

    # Space Y: The unit circle S^1. Its radius is 1.
    # The intrinsic metric on the circle corresponds to its circumference, l.
    r = 1.0
    l = 2 * math.pi * r

    # Step 2: State the problem and identify the parameters.
    print("Problem: Find the Gromov-Hausdorff distance between the interval [0, 1] and the unit circle S^1.")
    print(f"The interval is defined as [0, 1], so its length L = {L}.")
    print(f"The unit circle has radius r = {r}, so its circumference l = 2 * pi * r = {l:.5f}.")
    print("-" * 50)

    # The Gromov-Hausdorff distance d_GH(I_L, C_l) between an interval of length L
    # and a circle of circumference l depends on the comparison between L and the
    # circle's diameter (l/pi).

    # Step 3: Check the condition to select the correct formula.
    circle_diameter = l / math.pi
    print(f"We compare the interval length L = {L} with the circle's diameter l/pi = {circle_diameter:.5f}.")

    # According to the known theorem, we have two cases.
    if L <= circle_diameter:
        print(f"Since L <= l/pi ({L} <= {circle_diameter:.5f}), the formula is:")
        print("d_GH = (1/2) * sqrt(L^2 + (l/2)^2)")
        print("-" * 50)

        # Step 4: Calculate the distance using the selected formula.
        l_half = l / 2.0
        distance_squared = L**2 + l_half**2
        distance = 0.5 * math.sqrt(distance_squared)

        # Step 5: Print the final equation with substituted values and the result.
        print("Calculation Steps:")
        print(f"d_GH = (1/2) * sqrt({L:.1f}^2 + ({l:.5f}/2)^2)")
        print(f"d_GH = (1/2) * sqrt({L**2:.1f} + ({l_half:.5f})^2)")
        print(f"d_GH = (1/2) * sqrt({L**2:.1f} + {l_half**2:.5f})")
        print(f"d_GH = (1/2) * sqrt({distance_squared:.5f})")
        print(f"d_GH = 0.5 * {math.sqrt(distance_squared):.5f}")
        print("\nFinal Result:")
        print(f"The Gromov-Hausdorff distance is approximately: {distance:.5f}")

    else:  # This case, L > circle_diameter, does not apply here but is included for completeness.
        print(f"Since L > l/pi ({L} > {circle_diameter:.5f}), the formula is:")
        print("d_GH = (L + l/pi) / 2")
        print("-" * 50)
        
        # Calculate the distance using the other formula.
        distance = (L + circle_diameter) / 2.0
        
        print("Calculation Steps:")
        print(f"d_GH = ({L:.1f} + {l:.5f}/pi) / 2")
        print(f"d_GH = ({L:.1f} + {circle_diameter:.5f}) / 2")
        print(f"d_GH = {L + circle_diameter:.5f} / 2")
        print("\nFinal Result:")
        print(f"The Gromov-Hausdorff distance is approximately: {distance:.5f}")

if __name__ == '__main__':
    calculate_gromov_hausdorff_distance()