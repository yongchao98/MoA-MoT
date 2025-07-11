import math

def solve_container_problem():
    """
    This function solves the container design optimization problem.
    It calculates the costs for both box and cylinder designs and
    determines which is cheaper. It then prints the breakdown of the
    final cost calculation for the optimal design.
    """

    # --- Problem Constants ---
    REQUIRED_ENERGY = 1000  # MJ
    BALL_ENERGY = 25  # MJ
    BALL_COST = 1000  # USD
    BALL_DIAMETER = 4.0  # cm
    MATERIAL_COST_PER_CM2 = 200  # USD
    PI = math.pi

    # --- Step 1: Calculate minimum number of balls ---
    min_balls_needed = math.ceil(REQUIRED_ENERGY / BALL_ENERGY)

    # --- Step 2: Box Design Analysis ---
    # To hold at least 40 balls (min_balls_needed) and minimize surface area,
    # we need an arrangement n_l * n_w * n_h >= 40 with n_l, n_w, n_h as close
    # as possible. The combination (5, 4, 2) gives exactly 40 and is more
    # cube-like than other factorizations like (10, 2, 2) or (8, 5, 1).
    box_config = (5, 4, 2)
    box_num_balls = box_config[0] * box_config[1] * box_config[2]
    
    box_L = box_config[0] * BALL_DIAMETER
    box_W = box_config[1] * BALL_DIAMETER
    box_H = box_config[2] * BALL_DIAMETER
    
    box_surface_area = 2 * (box_L * box_W + box_L * box_H + box_W * box_H)
    box_container_cost = box_surface_area * MATERIAL_COST_PER_CM2
    box_ball_cost = box_num_balls * BALL_COST
    box_total_cost = box_container_cost + box_ball_cost

    # --- Step 3: Cylinder Design Analysis ---
    # We test different layer configurations. A layer of 7 balls (1 central, 6 hexagonal)
    # is very efficient.
    # Radius needed for 7 balls of radius 2cm is 6cm.
    # Height is determined by the number of layers needed.
    balls_per_layer = 7
    num_layers = math.ceil(min_balls_needed / balls_per_layer)
    cyl_num_balls = num_layers * balls_per_layer

    cyl_R = 6.0  # cm (fits 7 balls of radius 2cm)
    cyl_H = num_layers * BALL_DIAMETER
    
    # Check if dimensions are multiples of 0.5 cm
    # R=6.0, H=24.0, all OK.
    
    cyl_surface_area = (2 * PI * cyl_R**2) + (2 * PI * cyl_R * cyl_H)
    cyl_container_cost = cyl_surface_area * MATERIAL_COST_PER_CM2
    cyl_ball_cost = cyl_num_balls * BALL_COST
    cyl_total_cost = cyl_container_cost + cyl_ball_cost
    
    # --- Step 4: Compare and output the best design's cost calculation ---
    if box_total_cost < cyl_total_cost:
        print("The optimal design is a BOX container.")
        print("\n--- Final Cost Calculation ---")
        print(f"1. Ball Arrangement: {box_config[0]} x {box_config[1]} x {box_config[2]}")
        print(f"2. Total Balls: {box_num_balls}")
        print(f"3. Container Dimensions: Length = {box_L} cm, Width = {box_W} cm, Height = {box_H} cm")
        print("\nFinal Equation: C = (Number of Balls * Cost per Ball) + (Surface Area * Material Cost)")
        print(f"C = ({box_num_balls} * {BALL_COST}) + ({box_surface_area:.2f} * {MATERIAL_COST_PER_CM2})")
        print(f"C = {box_ball_cost} + {box_container_cost:.2f}")
        final_cost = box_total_cost

    else:
        print("The optimal design is a CYLINDER container.")
        print("\n--- Final Cost Calculation ---")
        print(f"1. Ball Arrangement: {num_layers} layers of {balls_per_layer} balls each")
        print(f"2. Total Balls: {cyl_num_balls}")
        print(f"3. Container Dimensions: Radius = {cyl_R} cm, Height = {cyl_H} cm")
        print("\nFinal Equation: C = (Number of Balls * Cost per Ball) + (Surface Area * Material Cost)")
        print(f"C = ({cyl_num_balls} * {BALL_COST}) + ({cyl_surface_area:.2f} * {MATERIAL_COST_PER_CM2})")
        print(f"C = {cyl_ball_cost} + {cyl_container_cost:.2f}")
        final_cost = cyl_total_cost

    print(f"C = {final_cost:.2f} USD")
    print(f"\n<<<{final_cost:.2f}>>>")


solve_container_problem()