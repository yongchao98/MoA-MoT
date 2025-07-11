import math

def solve_pioneer_probe_design():
    """
    Calculates the minimum cost to design a container for energy balls,
    comparing a box design and a more efficient cylindrical design.
    """

    # --- Constants ---
    ENERGY_PER_BALL_MJ = 25
    REQUIRED_ENERGY_MJ = 1000
    COST_PER_BALL_USD = 1000
    COST_PER_CM2_USD = 200
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4

    # --- Step 1: Calculate costs for the optimal BOX design ---
    # To minimize ball cost, we use exactly 40 balls.
    # To minimize surface area for a fixed volume, dimensions should be as close as possible.
    # Factors of 40 are explored: (2, 4, 5) is the most cube-like combination.
    n_balls_box = 40
    cost_balls_box = n_balls_box * COST_PER_BALL_USD

    # Dimensions based on 2x4x5 ball arrangement
    L = 2 * BALL_DIAMETER_CM  # 8 cm
    W = 4 * BALL_DIAMETER_CM  # 16 cm
    H = 5 * BALL_DIAMETER_CM  # 20 cm

    # Surface area and cost for the box
    area_box = 2 * (L*W + W*H + H*L)
    cost_container_box = area_box * COST_PER_CM2_USD
    total_cost_box = cost_balls_box + cost_container_box

    # --- Step 2: Calculate costs for the optimal CYLINDER design ---
    # We analyze a more efficient packing: hexagonal layers.
    # A layer of 7 balls (1 central, 6 around) is very efficient.
    balls_per_layer_cyl = 7
    min_balls_needed = math.ceil(REQUIRED_ENERGY_MJ / ENERGY_PER_BALL_MJ)

    # Calculate layers and total balls needed for this packing
    n_layers_cyl = math.ceil(min_balls_needed / balls_per_layer_cyl)
    n_balls_cyl = n_layers_cyl * balls_per_layer_cyl
    cost_balls_cyl = n_balls_cyl * COST_PER_BALL_USD

    # Dimensions for the cylinder
    # A 7-ball hexagonal cluster fits in a radius of 2 * ball_radius
    R_cyl = 2 * BALL_RADIUS_CM # 4 cm
    # Height is number of layers * ball diameter
    H_cyl = n_layers_cyl * BALL_DIAMETER_CM # 6 * 4 = 24 cm

    # Surface area and cost for the cylinder
    area_cyl = (2 * math.pi * R_cyl**2) + (2 * math.pi * R_cyl * H_cyl)
    cost_container_cyl = area_cyl * COST_PER_CM2_USD
    total_cost_cyl = cost_balls_cyl + cost_container_cyl

    # --- Step 3: Compare costs and determine the final answer ---
    if total_cost_cyl < total_cost_box:
        final_cost = total_cost_cyl
        winning_design = "Cylinder"
        num_balls = n_balls_cyl
        cost_ball_total = cost_balls_cyl
        cost_container_total = cost_container_cyl
        area = area_cyl
        radius = R_cyl
        height = H_cyl
    else:
        final_cost = total_cost_box
        winning_design = "Box"
        num_balls = n_balls_box
        cost_ball_total = cost_balls_box
        cost_container_total = cost_container_box
        area = area_box
        length = L
        width = W
        height = H

    print(f"The optimal design is a {winning_design} container.")
    print("\n--- Cost Calculation Breakdown ---")
    print(f"Total Cost = Ball Cost + Container Material Cost")
    
    # Print equation details
    print(f"\n1. Ball Cost:")
    print(f"   Number of balls = {num_balls}")
    print(f"   Cost per ball = ${COST_PER_BALL_USD}")
    print(f"   Ball Cost = {num_balls} * {COST_PER_BALL_USD} = ${cost_ball_total:,.2f}")

    print(f"\n2. Container Material Cost:")
    print(f"   Material cost = ${COST_PER_CM2_USD}/cm^2")
    if winning_design == "Cylinder":
        print(f"   Cylinder Radius (R) = {radius} cm")
        print(f"   Cylinder Height (H) = {height} cm")
        print(f"   Surface Area = 2 * pi * R^2 + 2 * pi * R * H")
        print(f"   Surface Area = 2 * {math.pi:.4f} * {radius}^2 + 2 * {math.pi:.4f} * {radius} * {height} = {area:.2f} cm^2")
    else: # Box
        print(f"   Box Length (L) = {length} cm")
        print(f"   Box Width (W) = {width} cm")
        print(f"   Box Height (H) = {height} cm")
        print(f"   Surface Area = 2 * (L*W + W*H + H*L)")
        print(f"   Surface Area = 2 * ({length}*{width} + {width}*{height} + {height}*{length}) = {area:.2f} cm^2")

    print(f"   Container Cost = {COST_PER_CM2_USD} * {area:.2f} = ${cost_container_total:,.2f}")

    print("\n--- Final Equation ---")
    print(f"Total Cost = ({num_balls} * {COST_PER_BALL_USD}) + ({COST_PER_CM2_USD} * {area:.2f})")
    print(f"Total Cost = ${cost_ball_total:,.2f} + ${cost_container_total:,.2f} = ${final_cost:,.2f}")
    
    # The final answer format as requested by the user prompt
    print(f"\n<<<{final_cost:.2f}>>>")

solve_pioneer_probe_design()