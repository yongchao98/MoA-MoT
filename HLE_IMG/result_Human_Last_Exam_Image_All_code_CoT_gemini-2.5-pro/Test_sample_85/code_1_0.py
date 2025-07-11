import math

def solve_cone_distance():
    """
    Calculates the furthest distance on a cone's surface from a point on the base.
    This function prints a step-by-step derivation of the solution.
    """
    # The variable 'd' is treated symbolically in the print statements.
    
    print("### Step-by-Step Solution ###\n")

    print("Step 1: Determine the properties of the cone.")
    print(" - The base diameter is 'd'.")
    print(" - The base radius (r) is d / 2.")
    print(" - The slant height (L), the distance from the apex to point P, is 'd'.\n")

    print("Step 2: Unroll the cone's surface into a 2D sector.")
    print(" - The radius of the unrolled sector is the slant height, L = d.")
    print(" - The arc length of the sector is the cone's base circumference, C = 2 * pi * r = 2 * pi * (d/2) = pi * d.")
    print(" - The angle (theta) of the sector is given by the formula: Arc Length = Radius * theta.")
    print(" - theta = (pi * d) / d = pi radians (or 180 degrees).\n")

    print("Step 3: Analyze the unrolled shape.")
    print(" - A sector with an angle of pi radians is a semicircle.")
    print(" - So, the cone's surface unrolls into a semicircle of radius 'd'.")
    print(" - Let's place the apex at the origin (0, 0) and the point P on the arc at (0, d).")
    print(" - We need to find the point Q(x, y) in the semicircle that is furthest from P(0, d).\n")

    print("Step 4: Locate the furthest point Q.")
    print(" - The distance is maximized when Q is on the boundary of the semicircle.")
    print(" - By checking the distance formula, the furthest points are the corners of the semicircle's diameter, located at (d, 0) and (-d, 0).\n")

    print("Step 5: Calculate the final distance.")
    print(" - We use the 2D distance formula between P(0, d) and Q(d, 0).")
    print("   Distance^2 = (x2 - x1)^2 + (y2 - y1)^2")
    print("   Distance^2 = (d - 0)^2 + (0 - d)^2")
    print("   Distance^2 = d^2 + d^2")
    
    final_coefficient_squared = 2
    print(f"   Distance^2 = {final_coefficient_squared} * d^2")

    print(f"\n   Distance = sqrt({final_coefficient_squared} * d^2)")
    print(f"   The furthest distance is d * sqrt({final_coefficient_squared}).\n")
    print("The numbers in the final equation are 1 (for each d^2 term) which sum to 2 under the square root.")

if __name__ == '__main__':
    solve_cone_distance()