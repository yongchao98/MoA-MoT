import sympy
from sympy import pi, sqrt

def solve_cone_distance():
    """
    Solves for the furthest distance on a cone's surface from a point P on the base.
    The code follows the logical steps outlined above.
    """
    # Use a symbol 'd' for the diameter, ensuring it's positive.
    d = sympy.Symbol('d', positive=True)

    print("Step 1: Define the cone's properties based on the problem statement.")
    slant_height = d
    base_radius = d / 2
    print(f"The slant height 's' is d.")
    print(f"The base radius 'r' is d/2.")
    
    print("\nStep 2: Determine the properties of the unrolled cone surface.")
    base_circumference = 2 * pi * base_radius
    sector_radius = slant_height
    # The angle theta = (arc length) / (radius) = (base circumference) / (slant height)
    sector_angle = base_circumference / sector_radius
    print(f"The unrolled surface is a sector of a circle with radius {sector_radius} and angle {sector_angle} radians.")
    print("This is a semicircle.")

    print("\nStep 3: Find the maximum distance on the unrolled semicircle.")
    # Place the apex at the origin (0,0).
    # Place point P at the top of the semicircle arc.
    P_coords = (0, d)
    # The furthest point Q is at a corner of the semicircle's diameter.
    Q_coords = (d, 0)
    print(f"In the unrolled plane, we place P at {P_coords} and the furthest point Q is at {Q_coords}.")

    print("\nStep 4: Calculate the distance between P and Q.")
    # Calculate the squared distance using the distance formula.
    dist_sq = (Q_coords[0] - P_coords[0])**2 + (Q_coords[1] - P_coords[1])**2
    # Take the square root to find the final distance.
    distance = sympy.sqrt(dist_sq)

    print("The distance is calculated using the distance formula: sqrt((d-0)^2 + (0-d)^2).")
    print(f"The squared distance is: 2*d**2")
    
    print("\nThe final equation for the distance, with each number shown explicitly, is:")
    # The final equation is distance = 1 * d * sqrt(2)
    print(f"Distance = {1} * d * sqrt({2})")

    print(f"\nThe furthest distance one can travel from P is: {distance}")

solve_cone_distance()