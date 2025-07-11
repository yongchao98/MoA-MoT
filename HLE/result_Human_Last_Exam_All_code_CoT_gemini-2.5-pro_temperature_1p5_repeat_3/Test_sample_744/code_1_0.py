import math

def solve_manufacturing_problem():
    """
    Calculates the highest value that can be made from the material by manufacturing balls.
    """
    # Step 1: Define the dimensions and values
    tube_outer_side = 20  # cm
    tube_thickness = 4    # cm
    tube_length = 100     # cm (1m)
    ball_radius = 2       # cm
    ball_diameter = ball_radius * 2

    value_whole_ball = 3
    value_welded_ball = 2

    # Derived dimension: the side length of the inner hollow square
    tube_inner_side = tube_outer_side - 2 * tube_thickness

    # Step 2: Calculate the number of balls for the "Whole Ball" strategy
    # The space (bounding box) required for one whole ball is a 4x4x4 cm cube.

    # First, calculate how many 4x4 squares can fit in the tube's cross-section.
    # The cross-section is a 20x20 square with a 12x12 hole.
    # This shape can be broken down into four 12x4 rectangles and four 4x4 corners.
    # Balls in the four 12x4 rectangular sections: 4 * (12cm_side / 4cm_diameter) = 12
    # Balls in the four 4x4 corner sections: 4 * (4cm_side / 4cm_diameter) = 4
    balls_per_cross_section = (4 * (tube_inner_side // ball_diameter)) + 4

    # Next, calculate how many 4cm layers fit along the tube's length.
    layers_for_whole_balls = tube_length // ball_diameter
    
    # Calculate total whole balls and their value
    total_whole_balls = balls_per_cross_section * layers_for_whole_balls
    total_value_from_whole_balls = total_whole_balls * value_whole_ball

    print("--- Strategy 1: Manufacturing Whole Balls ---")
    print(f"Number of balls per cross-section layer: {balls_per_cross_section}")
    print(f"Number of layers ({ball_diameter}cm thick) along the tube: {tube_length} / {ball_diameter} = {layers_for_whole_balls}")
    print(f"Total whole balls produced: {balls_per_cross_section} * {layers_for_whole_balls} = {total_whole_balls}")
    print(f"Total value from whole balls: {total_whole_balls} balls * {value_whole_ball} value/ball = {total_value_from_whole_balls}\n")


    # Step 3: Calculate the number of balls for the "Welded Ball" strategy
    # The space (bounding box) for a half-ball is a 4x4x2 cm block.

    # The number of half-balls per cross-section is the same as for whole balls.
    # Calculate how many 2cm layers fit along the tube's length.
    layers_for_half_balls = tube_length // ball_radius
    
    # Calculate total half-balls, the resulting welded balls, and their value.
    total_half_balls = balls_per_cross_section * layers_for_half_balls
    total_welded_balls = total_half_balls // 2
    total_value_from_welded_balls = total_welded_balls * value_welded_ball
    
    print("--- Strategy 2: Manufacturing Welded Balls ---")
    print(f"Number of half-balls per cross-section layer: {balls_per_cross_section}")
    print(f"Number of layers ({ball_radius}cm thick) for half-balls: {tube_length} / {ball_radius} = {layers_for_half_balls}")
    print(f"Total half-balls produced: {balls_per_cross_section} * {layers_for_half_balls} = {total_half_balls}")
    print(f"Total welded balls produced: {total_half_balls} / 2 = {total_welded_balls}")
    print(f"Total value from welded balls: {total_welded_balls} balls * {value_welded_ball} value/ball = {total_value_from_welded_balls}\n")

    # Step 4: Compare values and find the maximum
    highest_value = max(total_value_from_whole_balls, total_value_from_welded_balls)

    print("--- Conclusion ---")
    print(f"The highest value that can be made is: {highest_value}")

solve_manufacturing_problem()