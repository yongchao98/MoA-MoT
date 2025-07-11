import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # 1. Define the parameters of the two metric spaces.
    # The interval [0, 1] has length l = 1.
    l = 1.0
    # The unit circle is defined as having radius r = 1.
    # Its circumference (total length) L is 2 * pi * r.
    r = 1.0
    L = 2 * math.pi * r

    print("This script calculates the Gromov-Hausdorff distance between the interval [0,1] and the unit circle.")
    print(f"The interval has length l = {l}.")
    print(f"The unit circle (radius r={r}) has circumference L = {L:.5f}.")
    print("-" * 30)

    # 2. The calculation is based on a known formula from geometry, which depends on a parameter k.
    # k = L / (pi * l)
    k = L / (math.pi * l)

    print("The distance is found using a formula dependent on a parameter k = L / (pi * l).")
    print(f"k = ({L:.5f}) / (pi * {l}) = {k:.2f}")

    # 3. For k >= 2, the distance formula is d_GH = (l/2) * sqrt(1 + (k-1)^2).
    # Since our calculated k = 2.0, this is the correct formula to use.
    print("\nSince k=2, we use the formula for k >= 2:")
    print("d_GH = (l/2) * sqrt(1 + (k-1)^2)")
    print("-" * 30)

    # 4. Substitute the values into the formula and print the result.
    print("Substituting the values into the formula:")
    l_val = l
    k_val = k
    
    # We will output each step of the calculation as requested.
    part1_str = f"({l_val}/2)"
    k_minus_1_str = f"({k_val:.1f} - 1)"
    sqrt_inner_str = f"1 + {k_minus_1_str}^2"
    
    print(f"Distance = {part1_str} * sqrt({sqrt_inner_str})")

    k_minus_1_val = k_val - 1
    sqrt_inner_val = 1 + k_minus_1_val**2
    sqrt_val = math.sqrt(sqrt_inner_val)
    distance = (l_val / 2) * sqrt_val

    print(f"Distance = {l_val/2} * sqrt(1 + {k_minus_1_val:.1f}^2)")
    print(f"Distance = {l_val/2} * sqrt({sqrt_inner_val})")
    print(f"Distance = {l_val/2} * {sqrt_val:.5f}")
    
    print("\n" + "-"*30)
    print("The exact answer is sqrt(2)/2.")
    print(f"The final calculated Gromov-Hausdorff distance is: {distance:.5f}")


if __name__ == '__main__':
    calculate_gromov_hausdorff_distance()