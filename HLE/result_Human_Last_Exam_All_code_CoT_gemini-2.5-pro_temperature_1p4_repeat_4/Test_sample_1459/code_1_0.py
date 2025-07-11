import math

def calculate_gh_distance_segment_circle():
    """
    Calculates the Gromov-Hausdorff distance between the interval [0,1]
    and the unit circle S^1.
    """
    # Parameters for the interval [0, a]
    a = 1.0

    # Parameters for the unit circle (radius R=1)
    R = 1.0
    circumference = 2 * math.pi * R

    # Diameters of the two spaces
    diam_interval = a
    diam_circle = math.pi * R # Half the circumference

    print("This script calculates the Gromov-Hausdorff distance between the interval [0, 1] and the unit circle.")
    print("-" * 80)
    print(f"Space X: Interval [0, 1]")
    print(f"Diameter of X = {diam_interval}")
    print("")
    print(f"Space Y: Unit Circle (Radius {R})")
    print(f"Circumference of Y = 2 * pi * {R} = {circumference:.5f}")
    print(f"Diameter of Y = pi * {R} = {diam_circle:.5f}")
    print("-" * 80)

    # The GH distance between a segment of length 'a' and a circle of circumference 'C'
    # is (C/2 - a)/2, provided C >= 2a.
    # C/2 is the diameter of the circle.
    
    if circumference >= 2 * a:
        print("The condition C >= 2a holds, so we can use the formula: d_GH = (diam(Y) - diam(X)) / 2.")
        
        gh_distance = (diam_circle - diam_interval) / 2.0
        
        print("\nFinal Equation:")
        print(f"d_GH = (Diameter_Circle - Diameter_Interval) / 2")
        print(f"d_GH = ({diam_circle:.5f} - {diam_interval}) / 2")
        print(f"d_GH = {diam_circle - diam_interval:.5f} / 2")
        print(f"d_GH = {gh_distance:.5f}")

        print("\nSymbolically, the result is (pi - 1) / 2.")
        
    else:
        # This case is not reached for the given problem
        print("A different formula would be required as the condition C >= 2a is not met.")

calculate_gh_distance_segment_circle()
