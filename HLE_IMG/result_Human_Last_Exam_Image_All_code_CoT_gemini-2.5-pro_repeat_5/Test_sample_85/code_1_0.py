import math

def solve_cone_distance():
    """
    This function explains and calculates the furthest distance on a specific cone's surface.
    """

    print("Step 1: Determine the geometry of the unrolled cone surface.")
    print("The cone's base diameter is d, and its slant height is also d.")
    print("When unrolled, the surface forms a sector of a circle with radius d.")
    print("The sector angle is (pi * d) / d = pi radians (a semicircle).")
    print("\nStep 2: Find the furthest distance on the semicircle.")
    print("The problem is equivalent to finding the point Q in a semicircle of radius d")
    print("that is furthest from a point P on its arc.")
    print("The maximum distance is found to be from a point on the arc to a point on the opposite end of the diameter.")
    print("\nStep 3: State the final distance formula.")
    print("The maximum squared distance is 2 * d^2.")
    
    # The final equation is Distance = d * sqrt(2)
    # The number in this equation is 2.
    coefficient = 2
    sqrt_of_coefficient = math.sqrt(coefficient)
    
    print("\nThe final equation for the furthest distance is:")
    print(f"Distance = d * sqrt({coefficient})")
    print(f"Which simplifies to approximately d * {sqrt_of_coefficient:.4f}")
    
solve_cone_distance()