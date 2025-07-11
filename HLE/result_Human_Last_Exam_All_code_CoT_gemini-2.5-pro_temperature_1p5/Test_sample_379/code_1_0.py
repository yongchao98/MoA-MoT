import math

def solve_cube_locus():
    """
    Calculates the length of a locus of points on a cube's surface and expresses it as a percentage.
    
    The problem defines a cube with side length r, a point P at the midpoint of an edge,
    and a curve C which is the locus of points on the surface at a distance r from P.
    This script calculates the total length of C and divides it by 2 * pi * r.
    """

    # We express lengths in terms of pi and r. The value of r will cancel out in the final ratio.
    # The explanation of how each component is derived is in the text above.

    # 1. Length of the curve on the 2 faces adjacent to P's edge.
    # Each face contributes an arc of length (pi * r / 3).
    len_adj_faces = 2 * (math.pi / 3)  # This is the length in units of r

    # 2. Length on the next 2 faces (those adjacent to the first pair).
    # Each of these faces contributes a length of 2 * (pi * r / 3).
    len_secondary_faces = 2 * (2 * math.pi / 3) # This is the length in units of r

    # 3. Length on the final 2 "far" faces.
    # The locus on these faces consists of single points, so the length is 0.
    len_far_faces = 0

    # 4. Calculate the total length L of the curve C by summing the parts.
    # Let's represent the components of the sum clearly.
    # The equation for the total length L is the sum of the lengths on all faces.
    # L = (Length on Face_adj_1) + (Length on Face_adj_2) + ...
    # L = (pi*r/3) + (pi*r/3) + (2*pi*r/3) + (2*pi*r/3) + 0 + 0
    # L = 6 * pi * r / 3 = 2 * pi * r
    total_length_over_r = len_adj_faces + len_secondary_faces + len_far_faces

    # 5. The problem asks to compute the ratio L / (2 * pi * r).
    # Let's print the final equation with its components.
    numerator_val = 2 # The coefficient of pi*r in the total length
    denominator_val = 2 # The coefficient of pi*r in the divisor
    
    print("The total length of the curve C is L = 2 * π * r.")
    print("The problem asks for the ratio of L to (2 * π * r).")
    print("\nFinal Equation:")
    print(f"    ({numerator_val} * π * r)")
    print( "--------------------" )
    print(f"    ({denominator_val} * π * r)")

    # Perform the final division.
    ratio = total_length_over_r / (2 * math.pi)
    print(f"\nThe result of the division is {ratio:.1f}")

    # 6. Convert the ratio to a whole number percentage.
    percentage = ratio * 100
    print(f"\nAs a whole number percentage, the answer is {int(percentage)}%.")

    return int(percentage)

# Execute the function to get the final answer
final_answer = solve_cube_locus()
print(f"\n<<<{final_answer}>>>")