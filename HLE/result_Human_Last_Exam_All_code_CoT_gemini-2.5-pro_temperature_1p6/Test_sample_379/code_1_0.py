import math

def solve_cube_locus_problem():
    """
    Calculates the length of a specific locus on a cube's surface and expresses it as a percentage.

    The problem involves a cube with side length r and a point P at the midpoint of an edge.
    The locus C is the set of points on the cube's surface at a surface distance r from P.
    The length of C is found by summing the lengths of 6 identical circular arcs.
    """

    # Let the side length of the cube be r. We can set r=1 as it will cancel out.
    r = 1.0

    # The locus of points C is composed of 6 identical circular arcs.
    # We calculate the length of one arc.
    # The arc is a part of a circle with radius r. The central angle subtended by one arc is pi/3 radians (or 60 degrees).
    arc_radius = r
    arc_angle_radians = math.pi / 3
    
    # Length of one arc = radius * angle
    single_arc_length = arc_radius * arc_angle_radians
    
    # The total locus C consists of 6 such arcs.
    num_arcs = 6
    total_length_C = num_arcs * single_arc_length
    
    # The problem asks to divide the length of C by 2*pi*r.
    divisor = 2 * math.pi * r
    
    # Calculate the ratio.
    ratio = total_length_C / divisor
    
    # Convert the ratio to a whole number percentage.
    percentage = int(round(ratio * 100))
    
    # Print the steps of the final calculation as requested.
    print(f"The side length of the cube, r = {r}")
    print(f"The locus C consists of {num_arcs} arcs, each of length (pi * r / 3).")
    print(f"Total length of C = {num_arcs} * (pi * {r} / 3) = {total_length_C:.4f}")
    
    print("\nFinal Calculation:")
    print(f"The problem requires dividing the length of C by (2 * pi * r).")
    print(f"Value of C = {total_length_C:.4f}")
    print(f"Value of (2 * pi * r) = (2 * {math.pi:.4f} * {r}) = {divisor:.4f}")
    
    # Show the final equation with the numbers plugged in
    print(f"Result = ({total_length_C:.4f} / {divisor:.4f}) * 100")
    
    # Print the final answer.
    print(f"\nThe resulting ratio as a whole number percentage is: {percentage}%")

solve_cube_locus_problem()