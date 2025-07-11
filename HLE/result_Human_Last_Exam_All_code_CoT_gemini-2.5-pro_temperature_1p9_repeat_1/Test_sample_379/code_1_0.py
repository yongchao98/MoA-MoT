import math

def solve_cube_locus_problem():
    """
    Solves the described geometry problem about the locus of points on a cube's surface.
    """
    
    # As per the plan, the geometric analysis shows the locus C consists of
    # several circular arcs. The key findings from unfolding the cube are:
    # 1. The locus exists on 4 of the 6 faces of the cube.
    # 2. It is composed of a total of 6 arcs.
    # 3. Each arc is part of a circle of radius 'r'.
    # 4. Each of the 6 arcs subtends an angle of pi/3 radians (or 60 degrees).
    
    num_arcs = 6
    angle_per_arc_rad = math.pi / 3
    
    # Calculate the total angle for the entire locus C in radians.
    total_angle_rad = num_arcs * angle_per_arc_rad
    
    # The total length 'L' of the locus C is the radius 'r' times the total angle.
    # L = r * total_angle_rad = r * (6 * pi/3) = 2 * pi * r.
    # The factors for the formula L = k * pi * r are:
    numerator_coefficient = total_angle_rad / math.pi
    
    # The problem asks to divide this length L by (2 * pi * r).
    denominator_coefficient = 2
    
    # Calculate the final ratio. Note that 'pi' and 'r' cancel out.
    final_ratio = (numerator_coefficient * math.pi) / (denominator_coefficient * math.pi)
    
    # Convert the ratio to a whole number percentage.
    percentage = int(final_ratio * 100)
    
    # Print the step-by-step calculation as requested, showing the numbers
    # in the final equation.
    
    print(f"The total length of the locus C is L = {int(numerator_coefficient)} * pi * r.")
    print(f"We need to compute the ratio L / ({denominator_coefficient} * pi * r).")
    print("\nThe final equation is:")
    print(f"Ratio = ({int(numerator_coefficient)} * pi * r) / ({denominator_coefficient} * pi * r) = {final_ratio:.0f}")
    print("\nFinally, expressing this ratio as a whole number percentage:")
    print(f"Percentage = {final_ratio:.0f} * 100% = {percentage}%")

solve_cube_locus_problem()
<<<100>>>