import math

def solve_cube_locus_problem():
    """
    This function solves the geometry problem about the locus on a cube
    and explains the steps.
    """
    # Let r be the side length of the cube.
    # The point P is the midpoint of an edge.
    # The curve C is the locus of points on the surface at a distance r from P.

    # The total length of the curve C can be found by summing the lengths of
    # its constituent arcs. A detailed analysis of the cube's unfolded surface
    # reveals that the total length of C is exactly 2 * pi * r.

    # This can be understood by considering two families of paths from P, each
    # contributing to a part of the locus C:

    # 1. Paths towards the vertices of P's edge:
    #    The distance from P to a vertex is r/2. The path continues for another r/2,
    #    forming quarter-circles of radius r/2 on the two adjacent faces at that vertex.
    #    Length for 2 vertices = 2 * (2 * (1/4 * 2 * math.pi * (r/2))) = pi * r.

    # 2. Paths towards the centers of the faces containing P's edge:
    #    The distance from P to a face center is r/2. The path continues for another
    #    r/2, forming a semi-circle of radius r/2 on that face.
    #    Length for 2 faces = 2 * (1/2 * 2 * math.pi * (r/2)) = pi * r.

    # The total length of C is the sum of these lengths.
    # Total Length = (pi * r) + (pi * r) = 2 * pi * r.
    total_length_of_c_per_r = 2 * math.pi
    
    print("The total length of the locus C is calculated to be 2 * pi * r.")

    # The problem asks us to divide this length by 2 * pi * r.
    divisor_per_r = 2 * math.pi

    # The variable 'r' cancels out, so we can calculate the ratio directly.
    ratio = total_length_of_c_per_r / divisor_per_r

    # We output the numbers in the final equation as requested.
    print("\nWe are asked to calculate the ratio: (Length of C) / (2 * pi * r).")
    print("The equation is:")
    # Using f-string to format and display the numbers involved
    print(f"({2} * {math.pi:.5f} * r) / ({2} * {math.pi:.5f} * r) = {int(ratio)}")

    # Finally, we express this ratio as a whole number percentage.
    percentage = int(round(ratio * 100))

    print(f"\nThe resulting ratio is {int(ratio)}.")
    print(f"Expressed as a whole number percentage, the final answer is {percentage}%.")

# Execute the function to see the solution.
solve_cube_locus_problem()