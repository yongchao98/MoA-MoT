import math

def solve_manufacturing_problem():
    """
    Calculates the highest value that can be made from a square hollow tube
    by manufacturing spherical balls.
    """
    # --- Input Parameters ---
    outer_side_cm = 20
    length_m = 1
    thickness_cm = 4
    ball_radius_cm = 2
    whole_ball_value = 3
    welded_ball_value = 2

    # --- Calculations ---
    # Convert all units to cm
    length_cm = length_m * 100
    ball_diameter_cm = ball_radius_cm * 2

    # The problem can be solved by determining how many balls can be physically cut.
    # A ball with a 4cm diameter requires a 4x4x4 cm cube of material to be cut from.
    # The tube wall is 4cm thick, which is a perfect fit for the ball's diameter.

    # 1. Calculate how many 4x4 cm squares fit into the tube's cross-section.
    # The cross-section is a 20x20 square with a hollow center.
    # Inner side = Outer side - 2 * thickness
    inner_side_cm = outer_side_cm - 2 * thickness_cm

    # Number of 4x4 squares that fit along the outer side
    outer_squares_per_side = outer_side_cm // ball_diameter_cm
    # Number of 4x4 squares that fit in the hollow area's side
    inner_squares_per_side = inner_side_cm // ball_diameter_cm

    # Total 4x4 squares in the full 20x20 cross-section
    total_area_in_squares = outer_squares_per_side * outer_squares_per_side
    # Total 4x4 squares in the hollow 12x12 cross-section
    hollow_area_in_squares = inner_squares_per_side * inner_squares_per_side

    # Number of 4x4 squares in the material's cross-section
    material_squares_in_cross_section = total_area_in_squares - hollow_area_in_squares

    # 2. Calculate how many 4cm layers can be cut along the tube's length.
    num_layers = length_cm // ball_diameter_cm

    # 3. Calculate total number of whole balls.
    # This is the number of 4x4x4 cm cubes that can be cut from the tube.
    total_possible_whole_balls = material_squares_in_cross_section * num_layers

    # 4. Determine the maximum value.
    # For each 4x4x4cm block of material, we have two choices:
    # a) Make 1 whole ball, for a value of 3.
    # b) Make 2 half-balls, weld them into 1 welded ball, for a value of 2.
    # To maximize the total value, we should always choose option (a).
    # Therefore, the optimal strategy is to make as many whole balls as possible.

    max_value = total_possible_whole_balls * whole_ball_value

    print("Step 1: Calculate the number of balls that can be cut.")
    print(f"The tube's cross-section can fit {material_squares_in_cross_section} columns of balls.")
    print(f"The tube's length allows for {num_layers} balls per column.")
    print(f"Total number of whole balls = {material_squares_in_cross_section} * {num_layers} = {total_possible_whole_balls}")
    print("\nStep 2: Calculate the maximum value.")
    print("To maximize value, every possible ball should be made as a 'whole ball' (value 3) instead of a 'welded ball' (value 2).")
    print(f"The highest value is {total_possible_whole_balls} * {whole_ball_value} = {max_value}")

solve_manufacturing_problem()
<<<C>>>