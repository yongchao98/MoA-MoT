import math

def calculate_max_value():
    """
    Calculates the maximum value that can be obtained by manufacturing balls
    from a square hollow tube.
    """
    # Step 1: Define the dimensions and values from the problem.
    outer_side = 20  # cm
    length = 100     # cm
    thickness = 4    # cm
    ball_radius = 2  # cm
    whole_ball_value = 3
    welded_ball_value = 2

    # Step 2: Determine the size of the product for packing purposes.
    # A ball with a 2cm radius has a 4cm diameter.
    # To manufacture it, we need a cube of material with a side length equal to the diameter.
    ball_diameter = ball_radius * 2

    # Step 3: Calculate the usable cross-sectional area of the material.
    # The tube is hollow, so we subtract the inner area from the outer area.
    inner_side = outer_side - 2 * thickness
    outer_area = outer_side ** 2
    inner_area = inner_side ** 2
    material_cross_section_area = outer_area - inner_area

    # Step 4: Calculate how many balls can fit in the cross-section and along the length.
    # The cross-section of the cube needed for one ball is diameter x diameter.
    ball_packing_area = ball_diameter ** 2

    # Number of 'rows' of balls that fit into the material's cross-section.
    # The dimensions are perfect multiples, so we can use integer division.
    num_rows_in_cross_section = material_cross_section_area // ball_packing_area

    # Number of balls that fit along the length of the tube.
    num_balls_per_row = length // ball_diameter

    # Step 5: Since the value of a whole ball (3) is greater than a welded ball (2)
    # and the material is thick enough to make whole balls everywhere, the optimal
    # strategy is to only make whole balls.
    total_whole_balls = num_rows_in_cross_section * num_balls_per_row
    max_value = total_whole_balls * whole_ball_value

    # Print the step-by-step logic and the final equation.
    print("Step-by-step calculation:")
    print(f"1. The inner side of the hollow tube is {outer_side}cm - 2 * {thickness}cm = {inner_side}cm.")
    print(f"2. The usable cross-sectional area is {outer_side}^2 - {inner_side}^2 = {material_cross_section_area} cm^2.")
    print(f"3. Each ball has a {ball_diameter}cm diameter and requires a {ball_diameter}cm x {ball_diameter}cm = {ball_packing_area} cm^2 cross-section to be made.")
    print(f"4. Number of ball 'rows' that can fit in the cross-section: {material_cross_section_area} / {ball_packing_area} = {num_rows_in_cross_section} rows.")
    print(f"5. Number of balls that can fit along the {length}cm length: {length} / {ball_diameter} = {num_balls_per_row} balls per row.")
    print(f"6. Total number of whole balls that can be manufactured: {num_rows_in_cross_section} * {num_balls_per_row} = {total_whole_balls}.")
    print("\nSince a whole ball (value {0}) is more valuable than a welded ball (value {1}), we maximize value by making only whole balls.".format(whole_ball_value, welded_ball_value))
    print("\nFinal Equation:")
    print(f"{total_whole_balls} balls * {whole_ball_value} value/ball = {max_value}")

calculate_max_value()