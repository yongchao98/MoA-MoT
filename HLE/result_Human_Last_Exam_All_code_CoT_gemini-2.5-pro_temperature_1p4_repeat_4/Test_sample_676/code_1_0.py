import math

def solve_container_problem():
    """
    Calculates the minimum cost to build a container for energy balls,
    considering both box and cylinder shapes.
    """
    # --- Constants ---
    BALL_DIAMETER = 4.0
    BALL_COST = 1000.0
    MATERIAL_COST = 200.0
    MIN_BALLS_REQUIRED = 40
    PRECISION = 0.5

    min_overall_cost = float('inf')
    best_config = {}

    # --- Search for the optimal BOX container ---
    # We loop through integer choices for the number of balls along each axis.
    # The search limit of 41 is a heuristic that covers the minimum requirement (1x1x40)
    # and many other more cube-like, efficient configurations.
    search_limit = 41 
    for nx in range(1, search_limit):
        for ny in range(nx, search_limit): # Start from nx to avoid duplicate shapes (e.g., 2x3x4 vs 3x2x4)
            for nz in range(ny, search_limit): # Start from ny for the same reason
                
                num_balls = nx * ny * nz
                if num_balls < MIN_BALLS_REQUIRED:
                    continue

                length = nx * BALL_DIAMETER
                width = ny * BALL_DIAMETER
                height = nz * BALL_DIAMETER

                surface_area = 2 * (length * width + length * height + width * height)
                
                cost_of_balls = num_balls * BALL_COST
                cost_of_container = surface_area * MATERIAL_COST
                total_cost = cost_of_balls + cost_of_container

                if total_cost < min_overall_cost:
                    min_overall_cost = total_cost
                    best_config = {
                        "type": "Box",
                        "total_cost": total_cost,
                        "num_balls": num_balls,
                        "surface_area": surface_area,
                        "layout": f"{nx}x{ny}x{nz}",
                        "dims": f"L={length}cm, W={width}cm, H={height}cm"
                    }

    # --- Search for the optimal CYLINDER container ---
    # We loop through integer choices for layers and grid layout within layers.
    for num_layers in range(1, search_limit):
        for nx_layer in range(1, search_limit):
             for ny_layer in range(nx_layer, search_limit):

                balls_per_layer = nx_layer * ny_layer
                if balls_per_layer == 0: continue
                
                num_balls = num_layers * balls_per_layer
                if num_balls < MIN_BALLS_REQUIRED:
                    continue

                # Cylinder Height
                height = num_layers * BALL_DIAMETER

                # Calculate the radius required to encircle the rectangular layer grid
                rect_half_width = nx_layer * BALL_DIAMETER / 2
                rect_half_height = ny_layer * BALL_DIAMETER / 2
                
                # Distance from the center of the rectangle to the center of a corner ball
                dist_to_corner_center = math.sqrt(rect_half_width**2 + rect_half_height**2)
                
                # Add the ball's radius to get the minimum radius for the container
                required_radius = dist_to_corner_center + (BALL_DIAMETER / 2)
                
                # Enforce manufacturing precision by rounding up to the nearest 0.5 cm
                actual_radius = math.ceil(required_radius / PRECISION) * PRECISION

                # Calculate surface area and cost
                surface_area = (2 * math.pi * actual_radius * height) + (2 * math.pi * actual_radius**2)
                cost_of_balls = num_balls * BALL_COST
                cost_of_container = surface_area * MATERIAL_COST
                total_cost = cost_of_balls + cost_of_container
                
                if total_cost < min_overall_cost:
                    min_overall_cost = total_cost
                    best_config = {
                        "type": "Cylinder",
                        "total_cost": total_cost,
                        "num_balls": num_balls,
                        "surface_area": surface_area,
                        "layout": f"{num_layers} layers of {nx_layer}x{ny_layer} balls",
                        "dims": f"R={actual_radius}cm, H={height}cm"
                    }

    # --- Print the Final Result ---
    if not best_config:
        print(0)
        print("<<<0>>>")
    else:
        final_total_cost = int(round(best_config['total_cost']))
        final_ball_cost = int(round(best_config['num_balls'] * BALL_COST))
        final_container_cost = int(round(best_config['surface_area'] * MATERIAL_COST))
        
        print("Final Equation:")
        print(f"Total Cost = (Number of Balls * Cost Per Ball) + (Container Surface Area * Cost Per cm^2)")
        print(f"Total Cost = ({best_config['num_balls']} * ${int(BALL_COST)}) + ({best_config['surface_area']:.2f} cm^2 * ${int(MATERIAL_COST)})")
        print(f"Total Cost = ${final_ball_cost} + ${final_container_cost} = ${final_total_cost}")
        print(f"<<<{final_total_cost}>>>")

solve_container_problem()
