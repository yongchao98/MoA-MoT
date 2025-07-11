import math

def solve_manufacturing_problem():
    """
    Calculates the maximum value of balls manufactured from a hollow tube.
    """
    # Material and product dimensions in cm
    tube_outer_width = 20
    tube_length = 100  # 1m
    tube_thickness = 4
    ball_radius = 2

    # Value of products
    value_whole_ball = 3
    value_welded_ball = 2 # This is a distractor as no material is left over

    # --- Calculations ---

    # 1. Calculate derived dimensions
    ball_diameter = ball_radius * 2
    tube_inner_width = tube_outer_width - 2 * tube_thickness

    # 2. Calculate the number of balls from the two larger walls
    # Dimensions: 20cm x 100cm, thickness 4cm
    balls_len_wise = tube_length // ball_diameter
    balls_large_wall_width_wise = tube_outer_width // ball_diameter
    total_balls_from_large_walls = 2 * (balls_len_wise * balls_large_wall_width_wise)

    # 3. Calculate the number of balls from the two smaller walls
    # Dimensions: 12cm x 100cm, thickness 4cm
    balls_small_wall_width_wise = tube_inner_width // ball_diameter
    total_balls_from_small_walls = 2 * (balls_len_wise * balls_small_wall_width_wise)

    # 4. Calculate total balls and maximum value
    total_whole_balls = total_balls_from_large_walls + total_balls_from_small_walls
    max_value = total_whole_balls * value_whole_ball

    # --- Output Results ---
    print("Step 1: Calculate the number of balls from the two larger walls (20cm x 100cm).")
    print(f"Number of balls = 2 * (length / diameter) * (width / diameter)")
    print(f"Number of balls = 2 * ({tube_length} / {ball_diameter}) * ({tube_outer_width} / {ball_diameter}) = {total_balls_from_large_walls}")
    print("-" * 20)
    print("Step 2: Calculate the number of balls from the two smaller walls (12cm x 100cm).")
    print(f"Number of balls = 2 * (length / diameter) * (width / diameter)")
    print(f"Number of balls = 2 * ({tube_length} / {ball_diameter}) * ({tube_inner_width} / {ball_diameter}) = {total_balls_from_small_walls}")
    print("-" * 20)
    print("Step 3: Calculate total balls and the maximum value.")
    print(f"Total whole balls = Balls from large walls + Balls from small walls")
    print(f"Total whole balls = {total_balls_from_large_walls} + {total_balls_from_small_walls} = {total_whole_balls}")
    print("\nSince all material is used making higher-value whole balls, the final calculation is:")
    print(f"Highest Value = Total whole balls * Value per whole ball")
    print(f"Highest Value = {total_whole_balls} * {value_whole_ball} = {max_value}")

solve_manufacturing_problem()
<<<C>>>