import math

def calculate_max_value():
    # Step 1: Define dimensions
    outer_side = 20  # cm
    length = 100  # cm, which is 1m
    thickness = 4  # cm
    ball_radius = 2  # cm
    ball_diameter = ball_radius * 2

    value_whole_ball = 3
    value_welded_ball = 2

    print("--- Material and Product Specifications ---")
    print(f"Tube outer side: {outer_side} cm")
    print(f"Tube length: {length} cm")
    print(f"Tube wall thickness: {thickness} cm")
    print(f"Ball diameter: {ball_diameter} cm\n")

    # Step 2: Calculate packing capacity
    # Along the length
    balls_per_length = length // ball_diameter
    print("--- Packing Calculation ---")
    print(f"Number of balls that can fit along the 100cm length: {length} / {ball_diameter} = {int(balls_per_length)}")

    # In the cross-section
    # The wall is 4cm thick, which perfectly fits the 4cm ball diameter.
    # We can calculate the number of 4cm columns around the perimeter.
    # Outer perimeter sections: 20cm / 4cm = 5 balls per side
    # Inner side dimension: 20cm - 2 * 4cm = 12cm
    # Inner perimeter sections: 12cm / 4cm = 3 balls per side
    # Total columns = 2 outer sides + 2 inner sides
    # Note: This counts columns, not individual ball slots in the 2D cross-section.
    # A better way is to calculate the number of 4x4 squares that fit in the cross-sectional area.
    # Top strip: 20cm / 4cm = 5 columns
    # Bottom strip: 20cm / 4cm = 5 columns
    # Side strips (remaining part): (20cm - 4cm - 4cm) / 4cm = 12cm / 4cm = 3 columns per side
    cols_top_bottom = outer_side // ball_diameter
    inner_side = outer_side - 2 * thickness
    cols_sides = inner_side // ball_diameter
    total_columns = cols_top_bottom + cols_top_bottom + cols_sides + cols_sides
    print(f"Number of ball columns that fit in the cross-section: (20/4) + (20/4) + (12/4) + (12/4) = {cols_top_bottom} + {cols_top_bottom} + {cols_sides} + {cols_sides} = {total_columns}\n")

    # Step 3: Calculate total number of whole balls
    total_whole_balls = total_columns * balls_per_length
    print("--- Total Production Capacity ---")
    print(f"Total number of whole balls that can be manufactured: {total_columns} columns * {int(balls_per_length)} balls/column = {total_whole_balls}\n")

    # Step 4: Determine optimal strategy and calculate max value
    print("--- Value Maximization ---")
    print(f"Value of a directly manufactured ball: {value_whole_ball}")
    print(f"Value of a welded ball (from two half-balls): {value_welded_ball}")
    print("To maximize value, we should produce the item with the highest value.")
    print("Since the material's dimensions allow for perfect packing of whole balls with no waste, the optimal strategy is to only make whole balls.\n")

    max_value = total_whole_balls * value_whole_ball
    print("--- Final Calculation ---")
    print(f"The highest possible value is calculated by multiplying the number of whole balls by their value:")
    print(f"Final Equation: {total_whole_balls} * {value_whole_ball} = {max_value}")

calculate_max_value()