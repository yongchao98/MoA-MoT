import math

def solve_cube_locus_problem():
    """
    Calculates the length of the specified locus on a cube's surface and expresses it
    as a percentage of 2*pi*r.
    """
    # We can use r=1.0 for simplicity as it will cancel out in the final ratio.
    r = 1.0

    # 1. By unfolding the cube, we find the locus C is composed of several circular arcs.
    # Due to symmetry, there are 6 such arcs of equal length.
    # - 1 arc on the first face containing P.
    # - 1 arc on the second face containing P.
    # - 2 arcs on the "left" adjacent face.
    # - 2 arcs on the "right" adjacent face.
    num_arcs = 6

    # 2. Calculate the length of a single arc.
    # On an unfolded plane, the arc is part of a circle with radius r. The arc is
    # defined by the boundaries of a single face. The subtended angle is pi/3 radians.
    angle_per_arc = math.pi / 3  # This is pi/6 - (-pi/6)
    single_arc_length = r * angle_per_arc

    # 3. Calculate the total length of the curve C.
    total_length = num_arcs * single_arc_length

    # 4. The problem asks to divide the total length by 2*pi*r.
    denominator = 2 * math.pi * r
    ratio = total_length / denominator

    # 5. Express the ratio as a whole number percentage.
    percentage = int(round(ratio * 100))

    # Output the steps of the calculation as requested.
    print("The problem is to find the length of a curve C on a cube of side r.")
    print("The curve C is the locus of points at a surface distance r from P, the midpoint of an edge.")
    print(f"The curve C is composed of {num_arcs} identical circular arcs.")
    print(f"The length of each arc is given by the formula: radius * angle.")
    print(f"In our case, each arc has a length of r * (pi/3).")
    print(f"The total length L of the curve is {num_arcs} * (r * pi / 3) = 2 * pi * r.")
    print(f"We are asked to compute the ratio L / (2 * pi * r).")
    print(f"Ratio = (2 * pi * r) / (2 * pi * r) = {ratio:.1f}")
    print(f"Expressed as a whole number percentage, the answer is {percentage}%.")

solve_cube_locus_problem()
print(f"\n<<<100>>>")