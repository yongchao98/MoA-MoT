import math

def solve_manufacturing_problem():
    """
    Calculates the maximum value of balls manufactured from a hollow tube.
    This solution is based on a geometric packing interpretation, as the
    tube's wall thickness matches the ball's diameter.
    """
    # 1. Define dimensions
    outer_side = 20  # cm
    thickness = 4  # cm
    length = 100  # cm (1 meter)
    ball_radius = 2  # cm

    # Values
    value_whole_ball = 3
    value_welded_ball = 2 # This is lower, so we avoid it if possible

    # 2. Calculate dimensions for packing
    inner_side = outer_side - 2 * thickness
    ball_diameter = ball_radius * 2

    # 3. Calculate how many balls fit in the cross-section
    # The cross-section is a frame. We find its area.
    tube_cross_section_area = outer_side**2 - inner_side**2
    
    # Each ball occupies a square footprint on the cross-section
    ball_footprint_area = ball_diameter**2
    
    # Number of columns of balls that fit
    num_columns = tube_cross_section_area // ball_footprint_area

    # 4. Calculate how many balls fit along the length
    num_balls_per_column = length // ball_diameter

    # 5. Calculate total balls
    # Since the length is a perfect multiple of the diameter, there are no leftovers.
    total_balls = num_columns * num_balls_per_column

    # 6. Determine optimal strategy and calculate max value
    # Since whole balls have a higher value (3 > 2), we maximize value
    # by making only whole balls.
    max_value = total_balls * value_whole_ball

    # Print the reasoning and the final equation with numbers
    print(f"First, we calculate the number of balls that can be geometrically packed into the tube.")
    print(f"The tube has a cross-sectional area of ({outer_side}^2 - {inner_side}^2) = {tube_cross_section_area} cm^2.")
    print(f"Each ball has a diameter of {ball_diameter} cm and occupies a {ball_diameter}x{ball_diameter} cm footprint ({ball_footprint_area} cm^2).")
    print(f"Number of ball 'columns' that fit in the cross-section = {tube_cross_section_area} / {ball_footprint_area} = {num_columns}.")
    print(f"The tube is {length} cm long, so we can fit {length} / {ball_diameter} = {num_balls_per_column} balls along its length.")
    print(f"Total number of balls = {num_columns} columns * {num_balls_per_column} balls/column = {total_balls} balls.")
    print("\nTo maximize value, we produce only whole balls (value {value_whole_ball}) instead of welded balls (value {value_welded_ball}).")
    print("The final calculation for the maximum value is:")
    print(f"{max_value} = {total_balls} * {value_whole_ball}")

solve_manufacturing_problem()
<<<C>>>