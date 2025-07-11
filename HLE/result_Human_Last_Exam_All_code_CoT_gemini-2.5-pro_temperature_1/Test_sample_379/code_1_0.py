import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus of points on a cube's surface and expresses it as a percentage.
    """
    # Let r be the side length of the cube. We can use r=1.0 for simplicity,
    # as it will cancel out in the final ratio.
    r = 1.0
    pi = math.pi

    # The problem asks for the length of the locus C of points at surface-distance r
    # from P, the midpoint of an edge. We find this length by summing the lengths
    # of the arcs that make up C on each of the cube's 6 faces.

    # Let P be the midpoint of the edge between the TOP and FRONT faces.

    # 1. Locus on the two faces adjacent to P's edge (TOP and FRONT).
    # By unfolding the TOP and FRONT faces, we find the locus on this combined
    # surface is a circular arc. The portion of the arc on the TOP face has a
    # length of (pi * r / 3).
    length_on_top = (pi * r) / 3
    # By symmetry, the length of the arc on the FRONT face is identical.
    length_on_front = (pi * r) / 3

    # 2. Locus on the faces adjacent to the TOP/FRONT pair (LEFT and RIGHT).
    # To reach the RIGHT face, there are two distinct shortest-path unfoldings:
    # - Path 1: Across the TOP face. This creates an arc of length (pi * r / 3).
    # - Path 2: Across the FRONT face. This creates a second arc of length (pi * r / 3).
    length_on_right = (pi * r) / 3 + (pi * r) / 3
    # By symmetry, the length on the LEFT face is the same.
    length_on_left = (pi * r) / 3 + (pi * r) / 3

    # 3. Locus on the faces opposite to P's edge (BOTTOM and BACK).
    # The shortest surface distance from P to any point on the BOTTOM or BACK face is r.
    # This minimum distance is achieved at only a single point on each face.
    # A single point has a length of 0.
    length_on_bottom = 0.0
    length_on_back = 0.0

    # 4. The total length of C is the sum of the lengths on all six faces.
    # The final equation for the total length is the sum of these components.
    # L = L_top + L_front + L_right + L_left + L_bottom + L_back
    print("The final equation for the total length L is a sum of the lengths of 6 arcs:")
    print(f"L = (pi*r/3) + (pi*r/3) + (2*pi*r/3) + (2*pi*r/3) + 0 + 0")
    print("The non-zero coefficients for pi*r in the sum are:")
    print(f"1/3 (from TOP), 1/3 (from FRONT), 2/3 (from RIGHT), 2/3 (from LEFT)")
    
    total_length_C = length_on_top + length_on_front + length_on_right + length_on_left + length_on_bottom + length_on_back
    
    print(f"\nSymbolically, Total Length C = (6/3)*pi*r = 2*pi*r")

    # 5. Divide the length of C by 2*pi*r as requested.
    ratio = total_length_C / (2 * pi * r)

    # 6. Convert the ratio to a whole number percentage.
    percentage = int(round(ratio * 100))

    print(f"The final answer as a whole number percentage is: {percentage}")

solve_cube_locus_problem()
<<<100>>>