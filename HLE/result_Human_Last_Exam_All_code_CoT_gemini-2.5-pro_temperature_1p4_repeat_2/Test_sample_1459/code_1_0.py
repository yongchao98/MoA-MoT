import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # Step 1: Define the parameters of the two metric spaces.
    # The interval is [0, 1], so its length L is 1.
    L = 1.0
    # The unit circle has a radius R of 1.
    # The intrinsic metric is the arc length, so the total length is the circumference C.
    R = 1.0
    C = 2 * math.pi * R

    print("This script calculates the Gromov-Hausdorff distance between the interval [0,1] and the unit circle.")
    print("-------------------------------------------------------------------------------------------------")
    print(f"Step 1: Define the properties of the spaces.")
    print(f"The length of the interval [0,1] is L = {L}.")
    print(f"The circumference of the unit circle (radius {R}) is C = 2 * pi * R = {C:.5f}.")
    print("\n")

    # Step 2: Use the formula for the distance between an interval and a circle.
    # Formula: d_GH = (1/2) * max(C/2, |L - C/2|)
    print("Step 2: Apply the formula for the Gromov-Hausdorff distance d_GH = (1/2) * max(C/2, |L - C/2|).")
    
    # Calculate the components of the formula
    c_div_2 = C / 2
    abs_diff = abs(L - c_div_2)
    
    print(f"First, we calculate the two terms inside the max() function:")
    print(f"Term 1: C/2 = {C:.5f} / 2 = {c_div_2:.5f}")
    print(f"Term 2: |L - C/2| = |{L} - {c_div_2:.5f}| = {abs_diff:.5f}")
    print("\n")

    # Step 3: Find the maximum of the two terms and compute the final result.
    max_val = max(c_div_2, abs_diff)
    
    print(f"Step 3: Determine the maximum value and calculate the final distance.")
    print(f"Next, we find the maximum of these two terms:")
    print(f"max({c_div_2:.5f}, {abs_diff:.5f}) = {max_val:.5f}")
    
    distance = 0.5 * max_val
    
    print(f"\nFinally, we compute the distance by multiplying by 1/2:")
    print(f"Final Equation: d_GH = 0.5 * {max_val:.5f} = {distance:.5f}")
    print("-------------------------------------------------------------------------------------------------")
    print(f"The Gromov-Hausdorff distance is exactly pi/2, which is approximately {distance:.5f}.")

if __name__ == '__main__':
    calculate_gromov_hausdorff_distance()
