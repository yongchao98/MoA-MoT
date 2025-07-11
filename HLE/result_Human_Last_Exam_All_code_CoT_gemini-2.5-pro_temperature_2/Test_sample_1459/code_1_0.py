import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle with the intrinsic metric.
    """

    # 1. Define the properties of the first space: the interval X = [0,1].
    # The metric is d(a, b) = |a - b|.
    interval_length = 1.0
    
    # The diameter of the interval is its length.
    diam_X = interval_length

    # 2. Define the properties of the second space: the unit circle Y = S^1.
    # "Unit circle" implies the radius is 1. The metric is the shortest arc length.
    circle_radius = 1.0
    
    # The circumference of the circle is 2 * pi * R.
    circumference_Y = 2 * math.pi * circle_radius
    # The diameter of the circle (with the intrinsic metric) is half its circumference.
    diam_Y = circumference_Y / 2.0

    # 3. The Gromov-Hausdorff distance between these two spaces is known to be
    # half the difference of their diameters.
    # d_GH(X, Y) = |diam_Y - diam_X| / 2
    gh_distance = (diam_Y - diam_X) / 2.0
    
    # 4. Output the final equation with all its numeric components, as requested.
    pi_val = math.pi
    one_val = 1.0
    two_val = 2.0
    
    print("This script calculates the Gromov-Hausdorff distance between the interval [0,1] and the unit circle.")
    print("-" * 40)
    print(f"The diameter of the interval [0, 1] is: {diam_X}")
    print(f"The diameter of the unit circle (radius {circle_radius}) is pi: {diam_Y:.6f}")
    print("-" * 40)
    print("The Gromov-Hausdorff distance is given by the formula:")
    print("d_GH = (Diameter of Circle - Diameter of Interval) / 2")
    print("\nFinal equation and result:")
    print(f"({pi_val:.6f} - {one_val}) / {two_val} = {gh_distance:.6f}")

calculate_gromov_hausdorff_distance()