import math

def solve_manufacturing_problem():
    """
    Calculates the maximum value that can be obtained by manufacturing balls from a hollow tube.
    """
    # Step 1: Define the dimensions
    outer_side = 20  # cm
    thickness = 4    # cm
    length = 100     # cm, which is 1m
    ball_radius = 2  # cm

    # Step 2: Calculate the space required for one ball
    ball_diameter = ball_radius * 2
    print(f"The material is a square hollow tube with outer side {outer_side}cm, thickness {thickness}cm, and length {length}cm.")
    print(f"We want to make balls with a radius of {ball_radius}cm, which means a diameter of {ball_diameter}cm.")
    print(f"To cut one whole ball, the machine requires a cubic block of material measuring {ball_diameter}x{ball_diameter}x{ball_diameter} cm.")
    print("-" * 20)

    # Step 3: Calculate the number of balls that can be made
    inner_side = outer_side - 2 * thickness
    material_cross_section_area = outer_side**2 - inner_side**2
    ball_block_cross_section_area = ball_diameter**2

    num_blocks_in_cross_section = material_cross_section_area / ball_block_cross_section_area
    num_blocks_along_length = length / ball_diameter

    total_blocks = num_blocks_in_cross_section * num_blocks_along_length

    print("First, let's calculate the number of cutting blocks that fit in a cross-section of the tube.")
    print(f"The inner side of the tube is {outer_side}cm - 2 * {thickness}cm = {inner_side}cm.")
    print(f"The cross-sectional area of the material is {outer_side}^2 - {inner_side}^2 = {material_cross_section_area} cm^2.")
    print(f"The cross-sectional area of a cutting block is {ball_diameter}^2 = {ball_block_cross_section_area} cm^2.")
    print(f"Number of blocks per cross-section = {material_cross_section_area} / {ball_block_cross_section_area} = {int(num_blocks_in_cross_section)}.")
    print("\nNext, let's see how many blocks fit along the length of the tube.")
    print(f"Number of blocks along the length = {length}cm / {ball_diameter}cm = {int(num_blocks_along_length)}.")
    print(f"\nTotal number of available 4x4x4cm blocks = {int(num_blocks_in_cross_section)} * {int(num_blocks_along_length)} = {int(total_blocks)}.")
    print("-" * 20)

    # Step 4: Determine the optimal manufacturing strategy
    value_whole_ball = 3
    value_welded_ball = 2
    print("For each 4x4x4cm block, we can either make:")
    print(f"  1. A whole ball, with a value of {value_whole_ball}.")
    print(f"  2. A welded ball (from two half-balls), with a value of {value_welded_ball}.")
    print(f"To maximize value, we should always choose to make a whole ball, as {value_whole_ball} > {value_welded_ball}.")
    print("-" * 20)

    # Step 5: Calculate the maximum total value
    max_value = total_blocks * value_whole_ball
    print("The highest possible value is obtained by making only whole balls.")
    print("Final Calculation:")
    print(f"{int(total_blocks)} (balls) * {value_whole_ball} (value per ball) = {int(max_value)}")

solve_manufacturing_problem()
<<<C>>>