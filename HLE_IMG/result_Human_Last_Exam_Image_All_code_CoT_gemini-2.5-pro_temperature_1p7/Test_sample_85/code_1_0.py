import math

def calculate_furthest_distance_on_cone():
    """
    This script calculates the furthest distance from a point P on a cone's 
    base circumference, traveling only on the surface.

    Problem givens:
    - Right cone
    - Base diameter = d
    - Slant height (distance from apex to base circumference) = d
    """

    # Step 1: Characterize the unrolled cone surface.
    # The lateral surface of a cone unrolls into a sector of a circle.
    # The radius of this sector is the cone's slant height, L. Given L = d.
    # The arc length of this sector is the cone's base circumference, C.
    # The base radius is r = d / 2, so C = 2 * pi * r = 2 * pi * (d / 2) = pi * d.
    #
    # The angle (theta) of the sector is given by the formula: theta = Arc Length / Radius.
    # theta = (pi * d) / d = pi radians (or 180 degrees).
    # A sector with an angle of 180 degrees is a semicircle.
    # So, the unrolled surface of the cone is a semicircle of radius d.
    
    # Step 2: Find the furthest point on the unrolled surface (the semicircle).
    # Let's place the semicircle in a 2D coordinate system. The apex of the cone
    # is at the origin (0,0), and the diameter of the semicircle lies on the x-axis.
    # The starting point P is on the base of the cone, which corresponds to the arc of the semicircle.
    # Due to symmetry, we can place P at the 'top' of the arc, so P is at coordinates (0, d).
    #
    # We need to find the point Q within the semicircle that is furthest from P. This point
    # must lie on the boundary of the semicircle. The boundary consists of the arc and the diameter.
    # By inspection (or calculus), the furthest points from P=(0,d) are the endpoints of the
    # semicircle's diameter: Q1=(-d, 0) and Q2=(d, 0).

    # Step 3: Calculate the distance.
    # We calculate the straight-line distance between P=(0,d) and, for instance, Q=(d,0).
    # This distance D is the maximum possible geodesic distance from P on the cone.
    # Using the distance formula: D = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    
    # D^2 = (d - 0)^2 + (0 - d)^2
    # D^2 = d^2 + (-d)^2
    # D^2 = d^2 + d^2
    # D^2 = 2 * d^2
    # D = sqrt(2 * d^2) = d * sqrt(2)

    number_in_final_equation = 2
    
    print("The final calculation for the furthest distance D is as follows:")
    print("Based on the geometry, the squared distance D^2 is derived from the coordinates of the start point P=(0,d) and furthest point Q=(d,0) on the unrolled surface:")
    print("D^2 = (d - 0)^2 + (0 - d)^2")
    print("D^2 = d^2 + d^2")
    print(f"D^2 = {number_in_final_equation} * d^2")
    print(f"Taking the square root gives D = d * sqrt({number_in_final_equation}).")
    
calculate_furthest_distance_on_cone()