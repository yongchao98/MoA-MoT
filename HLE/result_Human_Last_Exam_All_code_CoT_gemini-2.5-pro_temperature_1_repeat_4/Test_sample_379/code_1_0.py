import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus on a cube's surface and expresses it as a percentage.
    """
    # Let the side length of the cube be r. We can set r=1 as it will cancel out.
    r = 1.0
    
    # The locus C consists of four circular arcs.
    # The radius of the circle for each arc is r.
    # The angle subtended by each arc is pi/3 radians (60 degrees).
    num_arcs = 4
    angle_rad = math.pi / 3
    
    # Length of a single arc
    arc_length = r * angle_rad
    
    # Total length of the locus C
    total_length_C = num_arcs * arc_length
    
    # The value to divide by, which is the circumference of a circle of radius r
    divisor = 2 * math.pi * r
    
    # Calculate the ratio
    ratio = total_length_C / divisor
    
    # Convert the ratio to a percentage and round to the nearest whole number
    percentage = round(ratio * 100)
    
    print("Step 1: The locus C is composed of 4 circular arcs.")
    print(f"Step 2: Each arc has a length calculated as radius * angle. Here, radius = {r} and angle = pi/3.")
    print(f"Step 3: The total length of the locus C is {num_arcs} * {arc_length:.4f} = {total_length_C:.4f}.")
    print(f"Step 4: The problem requires dividing this length by 2 * pi * r, which is {divisor:.4f}.")
    final_equation = f"({total_length_C:.4f} / {divisor:.4f}) * 100"
    print(f"Step 5: The final calculation is {final_equation} = {ratio*100:.2f}%.")
    print(f"Step 6: Rounded to the nearest whole number, the result is {percentage}%.")

solve_cube_locus_problem()
print("\n<<<67>>>")