import math

def solve_cube_locus_problem():
    """
    Calculates the length of a specific locus on a cube's surface and expresses it as a ratio.

    The problem involves finding the length of the locus of points C at a surface
    distance r from a point P, where P is the midpoint of an edge of a cube with side r.
    The derivation shows that the total length of C is composed of 6 circular arcs.

    Each arc has a radius of r and subtends an angle of pi/3 radians.
    """
    
    # We can assume r=1 without loss of generality, as it will cancel out.
    # However, to show the equation clearly, we will keep it as a variable.
    r_symbol = 'r'
    
    # The length of a single arc segment
    # Length = radius * angle
    # The radius of the arc is r and the angle is pi/3
    num_arc_segments = 6
    arc_angle_rad = "pi/3"
    single_arc_length_str = f"r * {arc_angle_rad}"
    
    # Total length of the curve C is the sum of the lengths of all segments
    total_length_c_str = f"{num_arc_segments} * ({single_arc_length_str}) = 2 * pi * r"
    
    # The problem asks to divide this length by 2*pi*r
    denominator_str = f"2 * pi * r"
    
    # Performing the calculation. For numerical computation, we can set r=1.
    r_val = 1.0
    total_length_c_val = 2 * math.pi * r_val
    denominator_val = 2 * math.pi * r_val
    ratio = total_length_c_val / denominator_val
    
    # Convert the ratio to a whole number percentage
    percentage = int(ratio * 100)
    
    print("Step 1: Calculate the total length of the curve C.")
    print(f"The curve C consists of {num_arc_segments} identical arcs.")
    print(f"Length of one arc = radius * angle = {r_symbol} * ({arc_angle_rad})")
    print(f"Total length of C = {total_length_c_str}")
    print("-" * 20)
    
    print("Step 2: Divide the length of C by 2*pi*r.")
    print(f"Ratio = (Total Length of C) / ({denominator_str})")
    # Show the final equation with values, assuming r=1
    print(f"Final Equation: ({2:.1f}*pi*{r_val:.1f}) / ({2:.1f}*pi*{r_val:.1f}) = {ratio:.1f}")
    print("-" * 20)
    
    print("Step 3: Convert the ratio to a whole number percentage.")
    print(f"Percentage = {ratio:.1f} * 100 = {percentage}%")

solve_cube_locus_problem()

print(f"\n<<<{int(100.0)}>>>")