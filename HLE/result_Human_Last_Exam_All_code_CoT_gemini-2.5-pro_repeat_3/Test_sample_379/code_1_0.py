import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus of points on a cube's surface and expresses it
    as a percentage of 2*pi*r.
    """
    # Let r be the side length of the cube.
    # The value of r does not affect the final ratio, so we can use r=1 for calculation.
    r = 1.0

    # The locus of points C at a surface distance r from P (midpoint of an edge)
    # is composed of 6 identical circular arcs on the faces of the cube.
    # Through geometric analysis by unfolding the cube's faces, each arc is found
    # to be part of a circle of radius r, subtending an angle of pi/3 radians.

    # Number of identical arcs
    num_arcs = 6

    # Angle subtended by each arc (in radians)
    arc_angle_rad = math.pi / 3

    # Length of a single arc = radius * angle
    single_arc_length = r * arc_angle_rad

    # Total length of the curve C is the sum of the lengths of the 6 arcs
    total_length_C = num_arcs * single_arc_length

    # The reference length is the circumference of a circle with radius r
    reference_length = 2 * math.pi * r

    # The final ratio is the total length of C divided by the reference length
    ratio = total_length_C / reference_length

    # The answer is the ratio expressed as a whole number percentage
    final_percentage = int(round(ratio * 100))

    # --- Outputting the results step-by-step ---
    print("Step 1: Determine the properties of the curve C.")
    print(f"The curve C is composed of {num_arcs} identical arcs.")
    print(f"The radius for each arc is the cube's side length, r = {r:.1f}")
    print(f"The angle for each arc is pi/3 radians, which is {math.degrees(arc_angle_rad):.0f} degrees.")
    print("-" * 20)

    print("Step 2: Calculate the total length of curve C.")
    print(f"The length of a single arc is r * (pi/3).")
    print(f"Total length of C = {num_arcs} * (r * pi/3) = 2 * pi * r.")
    print(f"Numerically, Total Length = {total_length_C:.4f}")
    print("-" * 20)

    print("Step 3: Calculate the final percentage.")
    print(f"The reference length is 2 * pi * r = {reference_length:.4f}")
    print(f"The ratio is (Total Length of C) / (Reference Length) = {ratio:.2f}")
    print(f"Expressed as a percentage: {ratio * 100:.0f}%")
    print("-" * 20)
    
    print("Final Equation Calculation:")
    # Using the symbolic relationship Total Length = 2*pi*r
    # The equation is (2 * pi * r) / (2 * pi * r) * 100
    # To show the numbers as requested:
    print(f"( {total_length_C:.4f} / {reference_length:.4f} ) * 100 = {final_percentage}")

solve_cube_locus_problem()
<<<100>>>