import math

def solve_cone_distance_problem():
    """
    This function provides a step-by-step derivation for finding the furthest
    distance from a point P on the base of a specific cone.
    """
    print("--- Step 1: Unrolling the Cone ---")
    print("To find the furthest distance on the cone's surface, we unroll it into a flat sector of a circle.")
    print("The given parameters are:")
    print(" - Base diameter = d, which means the base radius r = d/2.")
    print(" - Slant height (distance from apex to P) L = d.")
    print("\nThe radius of the unrolled sector is equal to the slant height, so Sector Radius = d.")
    print("The arc length of the sector is equal to the circumference of the cone's base.")
    print("Base Circumference C = 2 * pi * r = 2 * pi * (d/2) = pi * d.")
    
    print("\n--- Step 2: Calculating the Sector's Angle ---")
    print("The arc length S of a sector is given by S = (Sector Radius) * theta, where theta is the angle in radians.")
    print("Here, S = C, so we have: pi * d = d * theta.")
    print(f"Solving for theta, we get theta = pi radians, which is {math.pi / math.pi * 180} degrees.")
    print("This means the unrolled surface of the cone is a perfect semicircle with radius d.")

    print("\n--- Step 3: Setting up the 2D Coordinate System ---")
    print("Let's place the apex of the cone at the origin (0, 0) of a 2D plane.")
    print("The semicircle lies in the upper half-plane. We imagine cutting the cone along the line from the apex to point P.")
    print("In the unrolled plane, point P is represented by the two endpoints of the semicircle's diameter.")
    print("So, the locations of P are P1 = (d, 0) and P2 = (-d, 0).")

    print("\n--- Step 4: Finding the Furthest Point Q ---")
    print("The true distance on the cone from any point Q to P is the minimum of the straight-line distances from Q to P1 and P2 in our 2D plane.")
    print("We want to find the point Q in the semicircle that maximizes this minimum distance.")
    print("The point that is 'furthest' from both P1 and P2 must lie on the line of symmetry between them, which is the perpendicular bisector.")
    print("The perpendicular bisector of the segment from (-d, 0) to (d, 0) is the y-axis (x=0).")
    print("So, our point Q has coordinates (0, y). To maximize the distance from P1 or P2, we need to maximize y.")
    print("Since Q must be inside or on the semicircle, the maximum possible value for y is d.")
    print("Therefore, the furthest point Q is at the coordinates (0, d).")
    
    print("\n--- Step 5: Calculating the Furthest Distance ---")
    print("The distance is the straight-line distance from P1(d, 0) to Q(0, d).")
    print("Using the distance formula: D = sqrt((x2 - x1)^2 + (y2 - y1)^2)")
    print("D = sqrt((0 - d)^2 + (d - 0)^2)")
    print("D = sqrt((-d)^2 + d^2)")
    
    # Explicitly showing the numbers in the final equation
    coeff1_sq = 1
    coeff2_sq = 1
    total_coeff = coeff1_sq + coeff2_sq
    
    print(f"D = sqrt({coeff1_sq}*d^2 + {coeff2_sq}*d^2)")
    print(f"D = sqrt({total_coeff}*d^2)")
    print(f"D = d * sqrt({total_coeff})")

    print("\nThus, the furthest distance one can travel from point P is d times the square root of 2.")

solve_cone_distance_problem()