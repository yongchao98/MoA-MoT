import math

def solve_manufacturing_problem():
    """
    Calculates the maximum value of balls that can be manufactured from a hollow tube.
    """
    # Material dimensions
    outer_dim_cm = 20
    length_cm = 100
    thickness_cm = 4

    # Product dimensions
    ball_radius_cm = 2
    ball_diameter_cm = ball_radius_cm * 2

    # Product values
    whole_ball_value = 3
    welded_ball_value = 2

    print("Step 1: Calculate the usable cross-sectional area of the tube.")
    inner_dim_cm = outer_dim_cm - 2 * thickness_cm
    outer_area = outer_dim_cm ** 2
    inner_area = inner_dim_cm ** 2
    material_cross_section_area = outer_area - inner_area
    print(f"The usable cross-sectional area is ({outer_dim_cm}x{outer_dim_cm}) - ({inner_dim_cm}x{inner_dim_cm}) = {outer_area} - {inner_area} = {material_cross_section_area} cm^2.\n")

    print("Step 2: Determine how many balls can be cut from the material.")
    # Each ball is cut from a cube of material with side length equal to its diameter
    space_per_ball_cross_section = ball_diameter_cm ** 2
    space_per_ball_length = ball_diameter_cm

    # Number of balls that fit across the material's cross-section
    balls_per_layer = material_cross_section_area // space_per_ball_cross_section
    print(f"The number of balls per layer is the material area divided by the ball's cross-section area: {material_cross_section_area} / {space_per_ball_cross_section} = {balls_per_layer}.")

    # Number of layers that fit along the tube's length
    num_layers = length_cm // space_per_ball_length
    print(f"The number of layers along the tube's length is: {length_cm} / {space_per_ball_length} = {num_layers}.\n")
    
    # Total number of whole balls
    total_whole_balls = balls_per_layer * num_layers
    
    print("Step 3: Calculate the maximum value.")
    # Since whole balls (value 3) are more valuable than welded balls (value 2) and the material
    # can be perfectly utilized to create whole balls, we should not produce any welded balls.
    max_value = total_whole_balls * whole_ball_value
    print("To maximize value, we should only produce whole balls.")
    print(f"The final equation for the total value is: {total_whole_balls} * {whole_ball_value} = {max_value}.\n")
    print(f"The highest value that can be made is {max_value}.")

solve_manufacturing_problem()
<<<C>>>