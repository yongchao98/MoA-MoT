import math

def solve_puzzle():
    """
    Calculates the maximum value that can be obtained from manufacturing balls from a hollow tube.
    """
    # Step 1: Define the dimensions and values from the problem description.
    outer_side = 20  # cm
    length = 100     # cm (1m)
    thickness = 4    # cm
    ball_radius = 2  # cm
    
    whole_ball_value = 3
    welded_ball_value = 2

    # Step 2: Determine the key dimensions for calculation.
    # Each ball has a 4cm diameter, so it can be manufactured from a 4x4x4cm cube.
    ball_unit_side = ball_radius * 2
    
    # The hollow tube's material can be viewed as four rectangular slabs.
    # The inner side dimension is the outer side minus the thickness on both sides.
    inner_side = outer_side - 2 * thickness

    print("Step 1: Analyzing the material dimensions.")
    print(f"The hollow tube has an outer side of {outer_side}cm, an inner side of {inner_side}cm, a thickness of {thickness}cm, and a length of {length}cm.")
    print(f"Each ball requires a {ball_unit_side}cm x {ball_unit_side}cm x {ball_unit_side}cm block of material.")
    print("-" * 30)

    # Step 3: Calculate how many balls can be packed into the material.
    # The material consists of two large slabs and two smaller slabs.
    
    # Large slabs (top and bottom of the tube)
    # Number of balls along each dimension of one large slab
    balls_in_len = length // ball_unit_side
    balls_in_width_large = outer_side // ball_unit_side
    balls_in_thickness = thickness // ball_unit_side
    
    balls_from_two_large_slabs = 2 * (balls_in_len * balls_in_width_large * balls_in_thickness)

    # Small slabs (left and right sides of the tube)
    # Number of balls along the width of one small slab
    balls_in_width_small = inner_side // ball_unit_side
    
    balls_from_two_small_slabs = 2 * (balls_in_len * balls_in_width_small * balls_in_thickness)
    
    # Total number of balls
    total_balls = balls_from_two_large_slabs + balls_from_two_small_slabs

    print("Step 2: Calculating the total number of balls that can be made.")
    print(f"Balls from the 2 large slabs ({outer_side}x{thickness}x{length}cm): 2 * ({balls_in_width_large} * {balls_in_thickness} * {balls_in_len}) = {balls_from_two_large_slabs}")
    print(f"Balls from the 2 small slabs ({inner_side}x{thickness}x{length}cm): 2 * ({balls_in_width_small} * {balls_in_thickness} * {balls_in_len}) = {balls_from_two_small_slabs}")
    print(f"Total possible balls: {balls_from_two_large_slabs} + {balls_from_two_small_slabs} = {total_balls}")
    print("-" * 30)

    # Step 4: Calculate the maximum value.
    # To maximize value, we must always choose to make whole balls (value 3) over welded balls (value 2).
    max_value = total_balls * whole_ball_value

    print("Step 3: Calculating the maximum value.")
    print("To maximize value, every ball produced will be a whole ball.")
    print(f"The final calculation is:")
    print(f"{total_balls} balls * {whole_ball_value} value per ball = {max_value}")

solve_puzzle()
<<<C>>>