import math

def solve_max_value():
    """
    Calculates the maximum value that can be produced from a square hollow tube.
    """
    # Material dimensions in cm
    outer_side = 20
    length = 100  # 1m = 100cm
    thickness = 4

    # Product dimensions in cm
    ball_radius = 2
    ball_diameter = ball_radius * 2

    # Value of products
    value_whole_ball = 3
    value_welded_ball = 2

    # Step 1: Calculate the cross-sectional area of the material
    inner_side = outer_side - 2 * thickness
    cross_section_area = outer_side**2 - inner_side**2

    print(f"The material is a square hollow tube with an outer side of {outer_side}cm, a thickness of {thickness}cm, and a length of {length}cm.")
    print(f"The inner side is calculated as {outer_side} - 2 * {thickness} = {inner_side}cm.")
    print(f"The cross-sectional area of the material is {outer_side}^2 - {inner_side}^2 = {cross_section_area} cm^2.\n")

    # Step 2: Calculate how many balls can be cut
    # To cut a ball of diameter 4cm, we need a 4x4cm square profile along the length.
    ball_profile_area = ball_diameter**2
    
    # Since the thickness (4cm) and the ball diameter (4cm) are the same,
    # the material can be perfectly divided.
    num_balls_in_cross_section = cross_section_area / ball_profile_area

    # Calculate how many 4cm thick layers can be cut from the 100cm length
    num_layers = length / ball_diameter

    total_balls = num_balls_in_cross_section * num_layers

    print(f"Each ball has a diameter of {ball_diameter}cm and requires a {int(ball_diameter)}x{int(ball_diameter)}cm profile for cutting.")
    print(f"Number of balls that can fit in the cross-section: {int(cross_section_area)} / {int(ball_profile_area)} = {int(num_balls_in_cross_section)}.")
    print(f"Number of {int(ball_diameter)}cm layers that can be cut from the {length}cm length: {length} / {int(ball_diameter)} = {int(num_layers)}.\n")

    # Step 3: Calculate total possible balls
    print("The total number of balls that can be manufactured is:")
    print(f"(Balls per cross-section) * (Number of layers) = {int(num_balls_in_cross_section)} * {int(num_layers)} = {int(total_balls)} balls.\n")

    # Step 4: Maximize the value
    # Since the value of a whole ball (3) is greater than a welded ball (2),
    # we should produce only whole balls to maximize the total value.
    max_value = total_balls * value_whole_ball

    print("To maximize profit, we should only make whole balls, as their value ({value_whole_ball}) is higher than welded balls ({value_welded_ball}).")
    print("\nThe final calculation for the highest possible value is:")
    print(f"Total Value = (Total number of balls) * (Value per whole ball)")
    print(f"Total Value = {int(total_balls)} * {value_whole_ball} = {int(max_value)}")

solve_max_value()
<<<C>>>