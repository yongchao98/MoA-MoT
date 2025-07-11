import math

def solve_cube_locus_problem():
    """
    This function calculates the solution to the specified geometry problem.
    
    The problem asks for the length of a locus C on a cube's surface, divided by 2*pi*r,
    expressed as a whole number percentage.
    
    Let r be the side length of the cube.
    Let P be the midpoint of an edge.
    Let C be the locus of points at a surface distance r from P.
    
    The calculation proceeds as follows:
    1. The locus C is composed of circular arcs on 4 faces of the cube.
    2. Calculate the length of the arcs on each face in terms of pi and r.
    3. Sum the lengths to get the total length of C.
    4. Divide the total length by 2*pi*r.
    5. Convert the result to a whole number percentage.
    """
    
    # The length of the arc on each of the 2 faces adjacent to P's edge is (pi * r) / 3.
    # The total length of the two arcs on each of the 2 "corner" faces is (2 * pi * r) / 3.
    
    # We can represent the coefficients of (pi * r) for the arc lengths.
    coeff_adj_face = 1/3
    coeff_corner_face = 2/3
    
    # Total length of C is the sum of lengths on 2 adjacent and 2 corner faces.
    # Total C = 2 * (coeff_adj_face * pi * r) + 2 * (coeff_corner_face * pi * r)
    total_coeff = 2 * coeff_adj_face + 2 * coeff_corner_face
    
    # The final equation for the length of C is:
    # C = 2 * (1/3 * pi * r) + 2 * (2/3 * pi * r)
    # C = (2/3 + 4/3) * pi * r
    # C = (6/3) * pi * r
    # C = 2 * pi * r
    
    print("Step 1: Calculate the total length of the locus C.")
    print("The length of C is the sum of 4 arcs on the cube's surface.")
    print("Length C = 2 * (arc on adjacent face) + 2 * (arcs on corner face)")
    print(f"Length C = 2 * (1/3 * pi * r) + 2 * (2/3 * pi * r)")
    print(f"Length C = ({2 * coeff_adj_face:.2f} + {2 * coeff_corner_face:.2f}) * pi * r = {total_coeff:.2f} * pi * r")
    print("So, the total length of C is exactly 2 * pi * r.")
    
    # Now, divide the length of C by 2 * pi * r
    denominator_coeff = 2
    ratio = total_coeff / denominator_coeff
    
    print("\nStep 2: Divide the length of C by 2 * pi * r.")
    print(f"Ratio = (2 * pi * r) / (2 * pi * r)")
    print(f"Ratio = {int(ratio)}")
    
    # Finally, convert to a whole number percentage
    percentage = round(ratio * 100)
    
    print("\nStep 3: Express the ratio as a whole number percentage.")
    print(f"Percentage = {int(ratio)} * 100%")
    print(f"Final Answer = {percentage}%")

solve_cube_locus_problem()
<<<100>>>