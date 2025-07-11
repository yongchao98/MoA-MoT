import math

def solve_manufacturing_problem():
    """
    Calculates the maximum value of balls manufactured from a square hollow tube.
    """
    # Step 1: Define Material and Product Dimensions
    tube_outer_dim_cm = 20
    tube_length_cm = 100  # 1m
    tube_thickness_cm = 4
    ball_radius_cm = 2
    ball_diameter_cm = ball_radius_cm * 2

    # Define values
    whole_ball_value = 3
    welded_ball_value = 2

    print("--- Material and Product Specs ---")
    print(f"Tube Outer Dimension: {tube_outer_dim_cm} cm")
    print(f"Tube Length: {tube_length_cm} cm")
    print(f"Tube Wall Thickness: {tube_thickness_cm} cm")
    print(f"Ball Diameter: {ball_diameter_cm} cm")
    print("-" * 35 + "\n")

    # Step 2: Since whole ball value (3) > welded ball value (2),
    # we prioritize making whole balls.

    # Step 3: Calculate Packing Capacity
    
    # Calculate inner dimension of the tube
    tube_inner_dim_cm = tube_outer_dim_cm - 2 * tube_thickness_cm

    # Calculate the cross-sectional area of the tube material
    outer_area = tube_outer_dim_cm ** 2
    inner_area = tube_inner_dim_cm ** 2
    material_cross_section_area = outer_area - inner_area

    # Calculate how many balls fit in the cross-section.
    # We model this by packing 4x4cm squares.
    packing_square_area = ball_diameter_cm ** 2
    # Use integer division as we can't have partial balls in the cross-section grid
    num_balls_per_slice = material_cross_section_area // packing_square_area

    # Calculate how many slices fit along the length
    # The thickness of each slice corresponds to the ball diameter
    num_slices = tube_length_cm // ball_diameter_cm

    # Total number of whole balls
    total_whole_balls = num_balls_per_slice * num_slices

    print("--- Calculation Steps ---")
    print(f"1. Number of balls that fit in the tube's cross-section: {material_cross_section_area} cm^2 / {packing_square_area} cm^2 = {int(num_balls_per_slice)} balls")
    print(f"2. Number of ball 'slices' along the tube's length: {tube_length_cm} cm / {ball_diameter_cm} cm = {int(num_slices)} slices")
    print(f"3. Total number of whole balls: {int(num_balls_per_slice)} * {int(num_slices)} = {int(total_whole_balls)} balls")
    print("-" * 35 + "\n")
    
    # Step 4: Calculate Maximum Value
    max_value = total_whole_balls * whole_ball_value

    print("--- Final Value Calculation ---")
    print("Since all material can be used for whole balls, which have the highest value, the maximum value is:")
    # The final equation with each number printed out
    print(f"Final Equation: {int(total_whole_balls)} (balls) * {whole_ball_value} (value per ball) = {int(max_value)}")

solve_manufacturing_problem()