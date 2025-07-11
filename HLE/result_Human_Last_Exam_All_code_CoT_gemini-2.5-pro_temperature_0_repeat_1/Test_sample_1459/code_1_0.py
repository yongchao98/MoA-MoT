import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle S^1.
    """
    # Step 1: Define the parameters of the two metric spaces.
    # Space X: The interval [0, 1]. Its length is L.
    L = 1.0

    # Space Y: The unit circle S^1. Its circumference is C.
    # The unit circle has a radius r = 1.
    r = 1.0
    C = 2 * math.pi * r

    # Step 2: State the known formula for the Gromov-Hausdorff distance
    # between a line segment of length L and a circle of circumference C.
    # The formula is: d_GH = (1/2) * max(L, C/2)
    
    # Step 3: Calculate the components of the formula.
    # We need to compare the length of the interval L with half the
    # circumference of the circle, C/2.
    C_half = C / 2

    # Step 4: Determine the maximum of L and C/2.
    max_val = max(L, C_half)

    # Step 5: Calculate the final Gromov-Hausdorff distance.
    gh_distance = 0.5 * max_val

    # Step 6: Print the explanation and the final result, showing each number
    # in the equation.
    print("This script calculates the Gromov-Hausdorff distance between the interval [0, 1] and the unit circle S^1.")
    print("-" * 80)
    print("The distance is calculated using the formula: d_GH = (1/2) * max(L, C/2)")
    print("where L is the length of the interval and C is the circumference of the circle.")
    print("-" * 80)
    
    print(f"1. The length of the interval [0, 1] is L = {L:.4f}")
    
    print(f"2. The circumference of the unit circle is C = 2 * pi * {r} = {C:.4f}")
    
    print(f"3. Half the circumference is C/2 = {C_half:.4f}")
    
    print(f"4. We find the maximum of L and C/2: max({L:.4f}, {C_half:.4f}) = {max_val:.4f}")
    
    print(f"5. Finally, the Gromov-Hausdorff distance is (1/2) of this maximum value:")
    print(f"   d_GH = (1/2) * {max_val:.4f} = {gh_distance:.4f}")
    
    print("-" * 80)
    print(f"The exact value of the distance is pi / 2.")

if __name__ == '__main__':
    calculate_gromov_hausdorff_distance()