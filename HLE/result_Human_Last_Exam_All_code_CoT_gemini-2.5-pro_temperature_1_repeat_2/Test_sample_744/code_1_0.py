def solve_max_value():
    """
    Calculates the maximum value that can be obtained by manufacturing balls from a hollow tube.
    """
    # 1. Define the initial parameters from the problem description.
    outer_side = 20  # cm
    length = 100     # cm, which is 1m
    thickness = 4    # cm
    ball_radius = 2  # cm
    value_whole = 3
    value_welded = 2

    print("Step 1: Determine the dimensions of the material and the space required for one ball.")
    # Calculate dimensions needed for the packing calculation.
    ball_diameter = ball_radius * 2
    inner_side = outer_side - 2 * thickness
    print(f"The hollow tube has a length of {length} cm and a wall thickness of {thickness} cm.")
    print(f"Each ball has a diameter of {ball_diameter} cm, so it requires a {ball_diameter}x{ball_diameter}x{ball_diameter} cm space to be cut.\n")

    print("Step 2: Calculate how many balls can be packed into the material.")
    # Calculate how many balls fit along the length.
    balls_along_length = length // ball_diameter
    print(f"Along the tube's length of {length} cm, we can fit: {length} / {ball_diameter} = {balls_along_length} balls.")

    # Calculate how many rows of balls fit in the cross-section.
    material_cross_section_area = outer_side**2 - inner_side**2
    ball_cross_section_area = ball_diameter**2
    rows_in_cross_section = material_cross_section_area // ball_cross_section_area
    print(f"The material's cross-section area is ({outer_side}^2 - {inner_side}^2) = {material_cross_section_area} cm^2.")
    print(f"The cross-section area needed per ball is {ball_diameter}^2 = {ball_cross_section_area} cm^2.")
    print(f"Therefore, we can fit {material_cross_section_area} / {ball_cross_section_area} = {rows_in_cross_section} rows of balls in the cross-section.\n")

    # Calculate the total number of balls.
    total_whole_balls = rows_in_cross_section * balls_along_length
    print(f"The total number of whole balls that can be manufactured is {rows_in_cross_section} rows * {balls_along_length} balls per row = {total_whole_balls} balls.\n")

    print("Step 3: Determine the optimal strategy and calculate the maximum value.")
    # Since whole balls have a higher value and there's no leftover material, we only make whole balls.
    print(f"The value of a whole ball is {value_whole}, and a welded ball is {value_welded}.")
    print("Because making a whole ball provides more value, and the material can be cut perfectly without waste, we will only make whole balls.")
    
    # Calculate the final maximum value.
    max_value = total_whole_balls * value_whole
    print("\nThe final calculation for the maximum value is:")
    print(f"{total_whole_balls} * {value_whole} = {max_value}")

solve_max_value()