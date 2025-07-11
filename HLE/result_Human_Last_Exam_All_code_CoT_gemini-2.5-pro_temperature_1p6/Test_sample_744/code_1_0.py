def solve_max_value():
    """
    Calculates the maximum value of balls that can be manufactured from a hollow tube.
    """
    # --- Step 1: Define Parameters ---
    # Tube dimensions
    outer_side_cm = 20
    length_cm = 100  # 1m
    thickness_cm = 4

    # Ball dimensions
    ball_radius_cm = 2
    ball_diameter_cm = ball_radius_cm * 2

    # Values
    value_whole_ball = 3
    value_welded_ball = 2

    # --- Step 2: Determine Manufacturing Block Size ---
    # To cut a whole ball, we need a cube with a side equal to the ball's diameter.
    block_side_cm = ball_diameter_cm

    # --- Step 3: Calculate the number of blocks (balls) ---
    # Calculate the inner side of the hollow tube
    inner_side_cm = outer_side_cm - 2 * thickness_cm

    # Calculate the cross-sectional area of the material
    area_material_cross_section = outer_side_cm**2 - inner_side_cm**2

    # Calculate the cross-sectional area of one manufacturing block
    area_block_cross_section = block_side_cm**2

    # Calculate how many blocks can fit in one slice of the tube's cross-section
    num_blocks_per_layer = area_material_cross_section // area_block_cross_section

    # Calculate how many layers of blocks can be cut along the tube's length
    num_layers = length_cm // block_side_cm

    # The entire material can be used to create whole balls
    total_whole_balls = num_blocks_per_layer * num_layers

    # Since the dimensions divide perfectly, there's no leftover material for half-balls
    total_welded_balls = 0

    # --- Step 4: Calculate the Maximum Value ---
    # The optimal strategy is to make whole balls as they have a higher value.
    max_value = (total_whole_balls * value_whole_ball) + (total_welded_balls * value_welded_ball)
    
    # --- Final Output ---
    print(f"From a {length_cm} cm long tube with {outer_side_cm}x{outer_side_cm} cm outer dimensions and {thickness_cm} cm thick walls, we can manufacture balls of {ball_diameter_cm} cm diameter.")
    print(f"Total number of whole balls that can be made: {total_whole_balls}")
    print(f"Total number of welded balls that can be made: {total_welded_balls}")
    print("\nFinal Value Calculation:")
    print(f"({total_whole_balls} whole balls * {value_whole_ball} value/ball) + ({total_welded_balls} welded balls * {value_welded_ball} value/ball) = {max_value}")

solve_max_value()