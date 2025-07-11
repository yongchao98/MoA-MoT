import math

def solve_max_value():
    """
    Calculates the maximum value that can be obtained by manufacturing balls
    from a square hollow tube.
    """
    # --- Given Parameters ---
    outer_side = 20  # cm
    length = 100  # cm (1m)
    thickness = 4  # cm
    ball_radius = 2  # cm
    value_whole_ball = 3
    value_welded_ball = 2

    # --- Derived Dimensions ---
    ball_diameter = ball_radius * 2
    inner_side = outer_side - 2 * thickness

    # --- Step 1: Calculate how many 4x4cm "logs" can be cut ---
    # The material's cross-section is a frame. Since the ball diameter equals the
    # thickness, we can cut logs with a 4x4cm cross-section.
    material_cross_section_area = outer_side**2 - inner_side**2
    log_cross_section_area = ball_diameter**2
    total_logs = int(material_cross_section_area / log_cross_section_area)

    # --- Step 2: Calculate value from making only whole balls ---
    # Each whole ball requires a 4cm (diameter) length from a log.
    num_whole_balls_per_log = math.floor(length / ball_diameter)
    total_whole_balls = total_logs * num_whole_balls_per_log
    value_from_whole_balls = total_whole_balls * value_whole_ball

    # --- Step 3: Calculate value from making only welded balls ---
    # Each welded ball requires two half-balls. Each half-ball requires a 2cm (radius) length.
    num_half_balls_per_log = math.floor(length / ball_radius)
    total_half_balls = total_logs * num_half_balls_per_log
    total_welded_balls = math.floor(total_half_balls / 2)
    value_from_welded_balls = total_welded_balls * value_welded_ball

    # --- Step 4: Compare and output the best strategy's equation ---
    print("To find the highest value, we compare two strategies:\n")
    
    print("1. Manufacturing Whole Balls:")
    print(f"   - Number of balls: {total_whole_balls}")
    print(f"   - Value per ball: {value_whole_ball}")
    print(f"   - Total Value = {total_whole_balls} * {value_whole_ball} = {value_from_whole_balls}\n")

    print("2. Manufacturing Welded Balls:")
    print(f"   - Number of balls: {total_welded_balls}")
    print(f"   - Value per ball: {value_welded_ball}")
    print(f"   - Total Value = {total_welded_balls} * {value_welded_ball} = {value_from_welded_balls}\n")

    if value_from_whole_balls > value_from_welded_balls:
        max_value = value_from_whole_balls
        print("Conclusion: The highest value is obtained by manufacturing whole balls.")
        print(f"Highest Value = {total_whole_balls} * {value_whole_ball} = {max_value}")
    else:
        max_value = value_from_welded_balls
        print("Conclusion: The highest value is obtained by manufacturing welded balls.")
        print(f"Highest Value = {total_welded_balls} * {value_welded_ball} = {max_value}")

# Run the calculation and print the result
solve_max_value()