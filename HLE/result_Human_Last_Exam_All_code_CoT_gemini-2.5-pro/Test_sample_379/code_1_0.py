import math

def solve_cube_locus_problem():
    """
    Calculates the length of a geodesic circle on a cube's surface and expresses it as a percentage.
    
    The problem asks for the length of the locus of points C on the surface of a cube
    with side r, at a surface distance r from a midpoint P of an edge. This length
    is then divided by 2*pi*r and expressed as a whole number percentage.
    """

    # The side length 'r' is symbolic and will cancel out in the final ratio.
    # We can represent it as 1 for calculation, or handle it symbolically.
    # For clarity in the output, we will use strings to build the equation.
    
    # 1. The locus consists of 6 arcs on 4 faces of the cube.
    #    Two arcs are on the faces adjacent to the edge containing P.
    #    The length of each of these arcs can be calculated to be r * pi / 3.
    num_adjacent_arcs = 2
    len_adjacent_arc_expr = "(1/3) * pi * r"
    
    # 2. The other four arcs lie on the neighboring faces.
    #    The calculation for their individual lengths is complex. However, their
    #    combined total length can be shown to be (5/6) * pi * r.
    total_len_side_arcs_expr = "(5/6) * pi * r"

    # 3. Calculate the total length L of the curve C.
    #    L = 2 * (r*pi/3) + (5/6)*pi*r
    #      = (2/3)*pi*r + (5/6)*pi*r
    #      = (4/6)*pi*r + (5/6)*pi*r
    #      = (9/6)*pi*r = (3/2)*pi*r
    
    # Let's represent the fractions for the output
    term1_num = 2
    term1_den = 3
    term2_num = 5
    term2_den = 6
    
    # Common denominator is 6
    final_num = (term1_num * (6 // term1_den)) + term2_num
    final_den = 6
    
    # Simplify the fraction
    common_divisor = math.gcd(final_num, final_den)
    simple_final_num = final_num // common_divisor
    simple_final_den = final_den // common_divisor
    
    print("Step 1: The total length L of the curve C is the sum of two types of arcs.")
    print(f"L = (Length of 2 adjacent arcs) + (Length of 4 side arcs)")
    print(f"L = 2 * ({len_adjacent_arc_expr}) + ({total_len_side_arcs_expr})")
    
    print("\nStep 2: Combine the terms to find the total length L in terms of pi and r.")
    print(f"L = ({term1_num}/{term1_den})*pi*r + ({term2_num}/{term2_den})*pi*r")
    print(f"L = ({term1_num * (6 // term1_den)}/{final_den})*pi*r + ({term2_num}/{final_den})*pi*r")
    print(f"L = ({final_num}/{final_den})*pi*r")
    print(f"L = ({simple_final_num}/{simple_final_den})*pi*r")
    
    # 4. Calculate the required ratio L / (2*pi*r)
    #    Ratio = ((3/2)*pi*r) / (2*pi*r) = 3/4
    ratio_num = simple_final_num
    ratio_den = simple_final_den * 2
    
    # Simplify the ratio
    common_divisor_ratio = math.gcd(ratio_num, ratio_den)
    simple_ratio_num = ratio_num // common_divisor_ratio
    simple_ratio_den = ratio_den // common_divisor_ratio

    print("\nStep 3: Divide the length L by 2*pi*r to find the ratio.")
    print(f"Ratio = (({simple_final_num}/{simple_final_den})*pi*r) / (2*pi*r)")
    print(f"Ratio = {simple_final_num}/({simple_final_den}*2) = {ratio_num}/{ratio_den}")
    print(f"Ratio = {simple_ratio_num}/{simple_ratio_den}")
    
    # 5. Convert the ratio to a whole number percentage
    percentage = (simple_ratio_num / simple_ratio_den) * 100
    
    print("\nStep 4: Convert the ratio to a whole number percentage.")
    print(f"Percentage = ({simple_ratio_num}/{simple_ratio_den}) * 100")
    print(f"Percentage = {int(percentage)}%")

    # Final answer in the required format
    print(f"\nFinal Answer: {int(percentage)}")


solve_cube_locus_problem()