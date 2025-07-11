import math

def solve_manufacturing_value():
    """
    Calculates the highest value that can be made from a hollow tube by manufacturing balls.
    """
    # 1. Define dimensions
    tube_outer_side = 20  # cm
    tube_thickness = 4    # cm
    tube_length = 100     # cm (1m)
    ball_radius = 2       # cm
    
    whole_ball_value = 3
    welded_ball_value = 2

    ball_diameter = ball_radius * 2
    tube_inner_side = tube_outer_side - 2 * tube_thickness

    # 2. Strategy 1: Calculate value from whole balls
    # A whole ball is cut from a cube with side length equal to the ball's diameter.
    whole_ball_block_side = ball_diameter
    
    # Calculate how many 4x4 cubes fit in the cross-section
    outer_grid_side = tube_outer_side // whole_ball_block_side
    inner_grid_side = tube_inner_side // whole_ball_block_side
    blocks_per_layer = outer_grid_side**2 - inner_grid_side**2
    
    # Calculate how many 4cm thick layers fit along the tube's length
    num_layers_whole_ball = tube_length // whole_ball_block_side
    
    # Total whole balls that can be made
    total_whole_balls = blocks_per_layer * num_layers_whole_ball
    value_from_whole_balls = total_whole_balls * whole_ball_value

    # 3. Strategy 2: Calculate value from welded balls
    # A half-ball is cut from a 4x4x2cm block.
    half_ball_block_length = ball_radius
    
    # The number of 4x4 blocks in the cross-section is the same
    # Calculate how many 2cm thick layers fit along the tube's length
    num_layers_half_ball = tube_length // half_ball_block_length
    
    # Total half-balls that can be made
    total_half_balls = blocks_per_layer * num_layers_half_ball
    
    # Two half-balls make one welded ball
    total_welded_balls = total_half_balls // 2
    value_from_welded_balls = total_welded_balls * welded_ball_value

    # 4. Compare and print the result
    if value_from_whole_balls >= value_from_welded_balls:
        print("The highest value is obtained by making whole balls.")
        print("Calculation:")
        print(f"{total_whole_balls} * {whole_ball_value} = {value_from_whole_balls}")
        highest_value = value_from_whole_balls
    else:
        print("The highest value is obtained by making welded balls.")
        print("Calculation:")
        print(f"{total_welded_balls} * {welded_ball_value} = {value_from_welded_balls}")
        highest_value = value_from_welded_balls

solve_manufacturing_value()
<<<C>>>