import math

def solve_manufacturing_problem():
    """
    Calculates the maximum value of balls that can be manufactured from a hollow tube.
    """
    # Parameters of the hollow tube
    outer_side = 20  # cm
    length = 100  # cm (1m)
    thickness = 4  # cm

    # Parameters of the ball
    ball_radius = 2  # cm
    ball_diameter = ball_radius * 2

    # Value of balls
    value_whole_ball = 3
    value_welded_ball = 2

    # --- Step 1: Analyze the geometry ---
    # The inner side of the hollow tube's cross-section is calculated.
    inner_side = outer_side - 2 * thickness

    print("This problem is a 3D packing problem. We need to determine how many balls can be cut from the material.")
    print(f"The tube has an outer side of {outer_side}cm, an inner side of {inner_side}cm, and a length of {length}cm.")
    print(f"Each ball has a diameter of {ball_diameter}cm.")
    print("\nSince the tube's thickness ({thickness}cm) is equal to the ball's diameter ({ball_diameter}cm), we can cut balls perfectly from the material.")
    print("We can model the hollow tube as four rectangular blocks:")

    # --- Step 2: Calculate how many balls fit in each section ---

    # Section 1 & 2: Two large blocks from the top and bottom faces (20cm x 4cm x 100cm)
    lanes_in_large_face = outer_side // ball_diameter
    balls_per_lane = length // ball_diameter
    balls_in_one_large_section = lanes_in_large_face * balls_per_lane

    print(f"\n1. For the two large blocks ({outer_side}cm x {thickness}cm x {length}cm):")
    print(f"   Balls along width = {outer_side} / {ball_diameter} = {lanes_in_large_face}")
    print(f"   Balls along length = {length} / {ball_diameter} = {balls_per_lane}")
    print(f"   Total balls per large block = {lanes_in_large_face} * {balls_per_lane} = {balls_in_one_large_section}")

    # Section 3 & 4: Two smaller blocks from the side faces (12cm x 4cm x 100cm)
    lanes_in_small_face = inner_side // ball_diameter
    balls_in_one_small_section = lanes_in_small_face * balls_per_lane

    print(f"\n2. For the two small blocks ({inner_side}cm x {thickness}cm x {length}cm):")
    print(f"   Balls along width = {inner_side} / {ball_diameter} = {lanes_in_small_face}")
    print(f"   Balls along length = {length} / {ball_diameter} = {balls_per_lane}")
    print(f"   Total balls per small block = {lanes_in_small_face} * {balls_per_lane} = {balls_in_one_small_section}")


    # --- Step 3: Calculate the total number of whole balls ---
    total_balls = (2 * balls_in_one_large_section) + (2 * balls_in_one_small_section)

    print("\nCalculating the total number of whole balls:")
    print(f"Total balls = (2 * {balls_in_one_large_section}) + (2 * {balls_in_one_small_section}) = {total_balls}")

    # --- Step 4: Calculate the maximum value ---
    # Since whole balls (value 3) are more valuable than welded ones (value 2),
    # and there's no waste material forcing us to make half-balls, we only make whole ones.
    max_value = total_balls * value_whole_ball

    print("\nTo get the highest value, we should only make whole balls (value 3 each).")
    print("\nFinal Equation:")
    print(f"{total_balls} * {value_whole_ball} = {max_value}")

solve_manufacturing_problem()