import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus on a cube's surface and expresses it as a percentage.
    
    The locus C is the set of points at a surface distance r from a midpoint P of an edge
    on a cube with side length r.
    
    The total length of the curve C is the sum of its parts:
    1. A semicircle of radius r on the two faces adjacent to P.
    2. Four identical circular arcs of radius r on the next set of adjacent faces.
    """
    
    # Let r be the side length of the cube. We can compute the length in terms of r.
    # The problem asks for a ratio, so r will cancel out. We can set r=1 for simplicity.
    r = 1.0

    # Part 1: The curve on the two faces adjacent to the edge containing P.
    # When unfolded, these form a semicircle of radius r.
    # The length of a semicircle is pi * r.
    semicircle_length = math.pi * r
    
    # Part 2: The curve on the four faces adjacent to the first two.
    # By unfolding the cube appropriately, we find that the curve on each of these 
    # four faces is a circular arc of radius r that subtends an angle of pi/3 radians (60 degrees).
    # The length of one such arc is r * angle_in_radians.
    num_smaller_arcs = 4
    arc_angle_rad = math.pi / 3
    single_arc_length = r * arc_angle_rad
    four_arcs_length = num_smaller_arcs * single_arc_length
    
    # The total length of the curve C is the sum of these parts.
    total_length_C = semicircle_length + four_arcs_length
    
    # The problem asks to divide this length by 2 * pi * r.
    divisor = 2 * math.pi * r
    ratio = total_length_C / divisor
    
    # Finally, express the answer as a whole number percentage.
    percentage = round(ratio * 100)
    
    # The problem asks to output the numbers in the final equation.
    # Let's show the breakdown of the calculation.
    # L = (pi*r) + 4 * (pi*r/3) = (3*pi*r + 4*pi*r)/3 = 7*pi*r/3
    # Ratio = (7*pi*r/3) / (2*pi*r) = 7/6
    
    # We can use variables to show the numbers clearly.
    # Multiples of pi*r for each part of the length
    semicircle_multiple = 1.0
    smaller_arc_multiple = 1.0/3.0
    
    print("Step-by-step calculation:")
    print(f"Length of semicircle part = {semicircle_multiple} * pi * r")
    print(f"Length of the four other arcs = {num_smaller_arcs} * {smaller_arc_multiple:.2f} * pi * r")
    print(f"Total length C = ({semicircle_multiple} + {num_smaller_arcs} * {smaller_arc_multiple:.2f}) * pi * r = {semicircle_multiple + num_smaller_arcs * smaller_arc_multiple:.2f} * pi * r = 7/3 * pi * r")
    print(f"The value to be calculated is (Length of C) / (2 * pi * r)")
    print(f"Ratio = (7/3 * pi * r) / (2 * pi * r) = 7/6")
    print(f"As a percentage: (7/6) * 100 = { (7/6)*100 :.2f}%")
    print(f"Rounded to the nearest whole number, the answer is: {percentage}%")

solve_cube_locus_problem()
<<<117>>>