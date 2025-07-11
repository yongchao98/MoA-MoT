import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus of points on a cube's surface and expresses
    it as a specified percentage.
    """
    # Let r be the side length of the cube. We can work with r=1 and the result will be scalable.
    r = 1.0

    # Step 1: Explain the setup.
    # We unfold the cube. Let P, the midpoint of an edge, be at the center of our
    # coordinate system for calculation. The locus of points on the surface at
    # distance r from P corresponds to a circle of radius r on the unfolded 2D net.
    # Due to the cube's geometry, this circle is split into 6 identical arcs
    # across 6 different faces of the cube.

    # Step 2: Calculate the length of a single arc.
    # Each arc is a segment of a circle of radius r. The angle (theta) subtended
    # by one of these arcs can be shown to be pi/3 radians (or 60 degrees).
    # The length of an arc is given by radius * angle.
    radius = r
    angle_per_arc_rad = math.pi / 3
    single_arc_length = radius * angle_per_arc_rad

    # Step 3: Calculate the total length of the locus C.
    num_arcs = 6
    total_length_C = num_arcs * single_arc_length

    # Step 4: The problem asks us to compute (C / (2 * pi * r)) as a percentage.
    # The final equation is: (6 * (pi/3) * r) / (2 * pi * r)
    divisor = 2 * math.pi * r
    final_ratio = total_length_C / divisor
    percentage = round(final_ratio * 100)

    # Step 5: Output the components of the final equation and the result.
    print("The length of the locus C is derived from several segments.")
    print(f"Number of arc segments: {num_arcs}")
    print(f"Length of each arc segment: (pi/3) * r")
    print("The total length C is the sum of these segments.")
    print(f"Equation for C = {num_arcs} * (pi/3) * r = {total_length_C/math.pi:.1f} * pi * r")
    print(f"The value we divide by is: 2 * pi * r")
    print("\nThe final ratio is (2 * pi * r) / (2 * pi * r) = 1.0")
    print(f"Expressed as a whole number percentage, the answer is: {percentage}%")
    
    # Returning the numerical answer for the final formatted output
    return int(percentage)

final_answer = solve_cube_locus_problem()
# <<<100>>>