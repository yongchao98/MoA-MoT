import math

def solve_cube_locus_problem():
    """
    Calculates the length of a specific locus on a cube's surface and expresses it as a percentage.
    
    The problem asks for (Length of C) / (2 * pi * r), where C is the locus of points
    on the surface of a cube of side r at a surface distance r from a midpoint of an edge, P.
    """

    # Part 1: Length of the locus on the two faces adjacent to the edge containing P.
    # When unfolded, these two faces form a 2r x r rectangle. The locus of points
    # at distance r from P (midpoint of the long side) is a semicircle of radius r.
    # Length of a semicircle = pi * r.
    # We can represent the symbolic parts of the equation with strings.
    len_adj_str = "pi * r"
    
    # Part 2: Length of the locus on the four faces adjacent to the first pair.
    # For each of these 4 faces, the locus is an arc of a circle with radius r.
    # The angle of each arc is pi/3 radians (60 degrees).
    # Length of one arc = radius * angle = r * (pi / 3).
    # Total length for the four arcs = 4 * (pi * r / 3).
    num_secondary_arcs = 4
    len_secondary_arc_str = "pi * r / 3"
    len_secondary_total_str = f"{num_secondary_arcs} * ({len_secondary_arc_str})"

    # Part 3: Total length of the locus C.
    # Total Length = (Length on adjacent faces) + (Length on secondary faces)
    # L = pi*r + 4*pi*r/3 = 3*pi*r/3 + 4*pi*r/3 = 7*pi*r/3
    total_len_numerator = 3 + 4
    total_len_denominator = 3
    total_len_str = f"{total_len_numerator} * pi * r / {total_len_denominator}"

    print("Step 1: The length of the locus on the two faces adjacent to P is a semicircle of radius r.")
    print(f"   Length_1 = {len_adj_str}\n")
    
    print("Step 2: The length of the locus on each of the four secondary adjacent faces is an arc of length r * (pi/3).")
    print(f"   Length_2 = {len_secondary_total_str}\n")

    print("Step 3: The total length of the locus C is the sum of these parts.")
    print(f"   Total Length = ({len_adj_str}) + ({len_secondary_total_str})")
    print(f"   Total Length = {total_len_str}\n")

    # Part 4: Divide the total length by 2*pi*r.
    # Ratio = (7 * pi * r / 3) / (2 * pi * r)
    # The 'pi * r' terms cancel out.
    # Ratio = (7/3) / 2 = 7/6
    final_ratio_num = total_len_numerator
    final_ratio_den = total_len_denominator * 2
    ratio = final_ratio_num / final_ratio_den

    print("Step 4: Divide the total length by 2*pi*r to find the ratio.")
    print(f"   Ratio = ({total_len_str}) / (2 * pi * r)")
    print(f"   Ratio = {final_ratio_num} / {final_ratio_den}\n")

    # Part 5: Convert the ratio to a whole number percentage.
    percentage = ratio * 100
    final_percentage = round(percentage)

    print("Step 5: Convert the ratio to a whole number percentage.")
    print(f"   Percentage = ({final_ratio_num} / {final_ratio_den}) * 100")
    print(f"   Percentage = {percentage:.2f}%")
    print(f"   Rounded to the nearest whole number, the answer is {final_percentage}%\n")
    
    return final_percentage

# Run the solver and print the final answer in the required format.
final_answer = solve_cube_locus_problem()
print(f"<<<{final_answer}>>>")
