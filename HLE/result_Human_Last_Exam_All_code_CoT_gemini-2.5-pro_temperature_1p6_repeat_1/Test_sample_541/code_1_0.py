import math

def calculate_hull_areas():
    """
    Calculates and prints the convex hull area for different cut shapes,
    each with a total length of 1 meter.
    """
    L = 1.0  # Total length of the cut in meters

    print(f"--- Comparing Convex Hull Areas for Cuts of Total Length {L}m ---\n")

    # Case A: Symmetry of an Equilateral Triangle (D3)
    # The optimal shape is a Y-shape (tripod) with 3 arms at 120 degrees.
    print("Shape: Y-cut (Symmetry of an Equilateral Triangle)")
    num_arms_Y = 3
    r_Y = L / num_arms_Y
    # The convex hull is an equilateral triangle. Its side length 's' can be found
    # using the law of cosines on the triangle formed by two arms.
    # s^2 = r^2 + r^2 - 2*r*r*cos(120) = 2r^2 - 2r^2*(-0.5) = 3r^2
    s_Y = r_Y * math.sqrt(3)
    area_Y = (math.sqrt(3) / 4) * s_Y**2
    print(f"The cut has 3 arms of length {r_Y:.4f} m meeting at 120 degrees.")
    print(f"The convex hull is an equilateral triangle with side length {s_Y:.4f} m.")
    print("The area calculation is: (sqrt(3) / 4) * side^2")
    print(f"Area = (sqrt(3) / 4) * ({s_Y:.4f})^2 = {area_Y:.6f} m^2")
    print("-" * 20)

    # Case D: Symmetry of a Square (D4)
    # A '+'-shaped cut is better than a square boundary.
    print("Shape: +-cut (Symmetry of a Square)")
    num_arms_sq = 4
    r_sq = L / num_arms_sq
    # The convex hull is a square with diagonal 2*r.
    # Area = (1/2) * d1 * d2 = (1/2) * (2r)^2 = 2r^2
    area_sq = 2 * r_sq**2
    print(f"The cut has 4 arms of length {r_sq:.4f} m meeting at 90 degrees.")
    print(f"The convex hull is a square with diagonal length 2*r = {2*r_sq:.4f} m.")
    print("The area calculation is: 2 * arm_length^2")
    print(f"Area = 2 * ({r_sq:.4f})^2 = {area_sq:.6f} m^2")
    print("-" * 20)
    
    # Case E: Symmetry of a Regular Hexagon (D6)
    # A 6-pronged star cut.
    print("Shape: 6-Pronged Star (Symmetry of a Regular Hexagon)")
    num_arms_hex = 6
    r_hex = L / num_arms_hex
    # The convex hull is a regular hexagon with side length s=r.
    # Area = (3 * sqrt(3) / 2) * s^2
    area_hex = (3 * math.sqrt(3) / 2) * r_hex**2
    print(f"The cut has 6 arms of length {r_hex:.4f} m meeting at 60 degrees.")
    print(f"The convex hull is a regular hexagon with side length {r_hex:.4f} m.")
    print("The area calculation is: (3*sqrt(3)/2) * side^2")
    print(f"Area = (3*sqrt(3)/2) * ({r_hex:.4f})^2 = {area_hex:.6f} m^2")
    print("-" * 20)

    # Case G: Symmetry of a Circle (O(2))
    # Cut is the circumference of a circle.
    print("Shape: Circular Cut (Symmetry of a Circle)")
    # Circumference = L = 2 * pi * r => r = L / (2*pi)
    r_circ = L / (2 * math.pi)
    # Area = pi * r^2
    area_circ = math.pi * r_circ**2
    print(f"The cut is a circle with circumference {L:.4f} m and radius {r_circ:.4f} m.")
    print("The area calculation is: pi * radius^2")
    print(f"Area = pi * ({r_circ:.4f})^2 = {area_circ:.6f} m^2")
    print("-" * 20)

    print("\nConclusion: The Y-shaped cut, which has the symmetry of an equilateral triangle, yields the largest convex hull area.")

if __name__ == '__main__':
    calculate_hull_areas()