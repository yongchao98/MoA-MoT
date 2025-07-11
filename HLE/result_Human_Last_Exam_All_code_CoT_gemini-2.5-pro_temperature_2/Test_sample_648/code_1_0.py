import math

def solve_design_problem():
    """
    Solves the container design optimization problem.
    """
    # --- Problem Constants ---
    MIN_ENERGY_MJ = 1000
    ENERGY_PER_BALL_MJ = 30
    COST_PER_BALL = 1000
    MAX_SURFACE_AREA = 1000
    COST_PER_AREA = 200
    BALL_DIAMETER = 4.0
    PRECISION_STEP = 0.5

    # Step 1: Calculate minimum number of balls needed
    min_balls_needed = math.ceil(MIN_ENERGY_MJ / ENERGY_PER_BALL_MJ)
    
    # Initialize tracking variables for the best design
    min_sa = float('inf')
    best_config = None

    # Set a practical search limit for ball grid dimensions
    search_limit = min_balls_needed + 5  # Search up to 39x39x39 grid

    # --- Search for the best BOX container ---
    for nx in range(1, search_limit):
        for ny in range(nx, search_limit):
            for nz in range(ny, search_limit):
                capacity = nx * ny * nz
                if capacity < min_balls_needed:
                    continue

                # Dimensions are multiples of ball diameter (4cm), which is a multiple of the 0.5cm precision
                L = nx * BALL_DIAMETER
                W = ny * BALL_DIAMETER
                H = nz * BALL_DIAMETER

                sa = 2 * (L * W + L * H + W * H)
                
                # If SA is already larger than a potential better solution or the max allowed, no need to continue with this nx, ny
                if sa > min_sa or sa > MAX_SURFACE_AREA:
                    break 
                
                min_sa = sa
                best_config = {
                    "type": "box",
                    "SA": sa,
                    "dims": (L, W, H),
                    "ball_grid": (nx, ny, nz),
                    "capacity": capacity,
                }
    
    # --- Search for the best CYLINDER container ---
    # In this search, nz is the number of layers, nx/ny define the layer arrangement
    for nx in range(1, search_limit):
        for ny in range(nx, search_limit):
            balls_per_layer = nx * ny
            if balls_per_layer == 0:
                continue
                
            nz = math.ceil(min_balls_needed / balls_per_layer)
            capacity = balls_per_layer * nz
            
            # Height is a multiple of ball diameter
            H = nz * BALL_DIAMETER
            
            # Radius must be large enough for the rectangular arrangement of balls in a layer
            # and quantized to the precision step
            min_internal_radius = 0.5 * math.sqrt((nx * BALL_DIAMETER)**2 + (ny * BALL_DIAMETER)**2)
            R = math.ceil(min_internal_radius / PRECISION_STEP) * PRECISION_STEP
            
            sa = 2 * math.pi * R * H + 2 * math.pi * R**2
            
            if sa <= MAX_SURFACE_AREA and sa < min_sa:
                min_sa = sa
                best_config = {
                    "type": "cylinder",
                    "SA": sa,
                    "dims": (R, H),
                    "ball_grid": (nx, ny, nz),
                    "capacity": capacity,
                }
                
    # --- Final Output ---
    if best_config:
        total_cost = best_config["SA"] * COST_PER_AREA + min_balls_needed * COST_PER_BALL
        sa_val = best_config["SA"]
        cost_area_val = COST_PER_AREA
        balls_val = min_balls_needed
        cost_ball_val = COST_PER_BALL
        
        print(f"A solution was found with a {best_config['type']} container.")
        if best_config['type'] == 'box':
            print(f"Optimal Configuration: Box with dimensions {best_config['dims']} cm for a ball grid of {best_config['ball_grid']}.")
        else:
            print(f"Optimal Configuration: Cylinder with Radius={best_config['dims'][0]} cm and Height={best_config['dims'][1]} cm for a ball grid of {best_config['ball_grid']}.")

        print(f"The container has a surface area of {sa_val:.2f} cm^2 and a capacity of {best_config['capacity']} balls.")
        print("To minimize cost, we will use the minimum required 34 balls.")
        print("\nThe final cost calculation is:")
        print(f"C = Surface Area * Cost per Area + Number of Balls * Cost per Ball")
        print(f"C = {sa_val:.2f} * {cost_area_val} + {balls_val} * {cost_ball_val} = {total_cost:.2f}")
        print(f"<<<{total_cost:.2f}>>>")

    else:
        # If no configuration meets the constraints
        print("No container design could be found that satisfies all the constraints (Surface Area <= 1000 cm^2 and Capacity >= 34 balls).")
        print(0)
        print("<<<0>>>")

solve_design_problem()