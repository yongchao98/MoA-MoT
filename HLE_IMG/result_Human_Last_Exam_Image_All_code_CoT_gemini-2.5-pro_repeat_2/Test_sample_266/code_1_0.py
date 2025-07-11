import math

def calculate_surface_area():
    """
    Calculates the surface area of the white shape based on a geometric model.
    The shape is modeled as a composite of 3 identical kites.
    """

    # The length of the first diagonal (d1) of the kite is the distance from the
    # center of the shape to an outer vertex, which corresponds to the side
    # length of the equilateral triangles in the background grid.
    d1 = 3.0

    # The second diagonal (d2) is calculated using the Law of Cosines on the
    # isosceles triangle formed by two edges of the white shape.
    # The length of these edges (s) is 2.5.
    s = 2.5
    # The angle (alpha) between these edges is 120 degrees.
    alpha_deg = 120.0
    alpha_rad = math.radians(alpha_deg)

    # Law of Cosines: d2^2 = s^2 + s^2 - 2*s*s*cos(alpha)
    d2_squared = s**2 + s**2 - 2 * s * s * math.cos(alpha_rad)
    d2 = math.sqrt(d2_squared)

    # The area of a single kite is 0.5 * d1 * d2
    area_kite = 0.5 * d1 * d2

    # The total area is the sum of the areas of the three kites.
    total_area = 3 * area_kite
    
    # Printing the equation with the final numbers
    print(f"Area = 3 * (0.5 * d1 * sqrt(s^2 + s^2 - 2*s*s*cos(alpha)))")
    print(f"Area = 3 * (0.5 * {d1} * sqrt({s}^2 + {s}^2 - 2*{s}*{s}*cos({alpha_deg})))")
    print(f"Area = 3 * (0.5 * {d1} * {d2:.4f})")
    print(f"Area = {total_area:.2f}")

calculate_surface_area()