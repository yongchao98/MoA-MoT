import math

def solve_container_problem():
    """
    Solves the container design optimization problem to find the minimum cost.
    """
    # --- Problem Constants ---
    ENERGY_REQ_MJ = 1000
    ENERGY_PER_BALL_MJ = 25
    COST_PER_BALL_USD = 1000
    BALL_DIAMETER_CM = 4
    MATERIAL_COST_PER_CM2_USD = 200

    # 1. Calculate minimum number of balls needed
    min_balls_needed = math.ceil(ENERGY_REQ_MJ / ENERGY_PER_BALL_MJ)

    min_total_cost = float('inf')
    best_config = None

    # 2. Search for the optimal box configuration
    # We check for total ball counts starting from the minimum required.
    # A slightly higher ball count might yield a more cubic shape,
    # potentially reducing surface area. We set a reasonable upper limit for the search.
    search_limit = min_balls_needed * 2
    for n_balls in range(min_balls_needed, search_limit):
        # Find all integer factor triplets (n_l, n_w, n_h) for n_balls
        # This represents the grid of balls inside the container
        for n_h in range(1, n_balls + 1):
            if n_balls % n_h == 0:
                area = n_balls // n_h
                for n_w in range(1, int(math.sqrt(area)) + 1):
                    if area % n_w == 0:
                        n_l = area // n_w

                        # We have a valid ball grid: n_l x n_w x n_h
                        # Calculate container dimensions
                        L = n_l * BALL_DIAMETER_CM
                        W = n_w * BALL_DIAMETER_CM
                        H = n_h * BALL_DIAMETER_CM

                        # Calculate surface area
                        surface_area = 2 * (L * W + W * H + H * L)

                        # Calculate total cost
                        cost_balls = n_balls * COST_PER_BALL_USD
                        cost_material = surface_area * MATERIAL_COST_PER_CM2_USD
                        total_cost = cost_balls + cost_material

                        # Check if this is the best solution so far
                        if total_cost < min_total_cost:
                            min_total_cost = total_cost
                            best_config = {
                                "n_balls": n_balls,
                                "grid": (n_l, n_w, n_h),
                                "dims": (L, W, H),
                                "surface_area": surface_area,
                                "cost_balls": cost_balls,
                                "cost_material": cost_material,
                                "total_cost": total_cost
                            }

    # 3. Output the results of the best configuration found
    if best_config:
        print("--- Optimal Design Found: Box Container ---")
        
        # Unpack for clarity
        n_balls = best_config["n_balls"]
        grid = best_config["grid"]
        dims = best_config["dims"]
        surface_area = best_config["surface_area"]
        cost_balls = best_config["cost_balls"]
        cost_material = best_config["cost_material"]
        total_cost = best_config["total_cost"]

        print(f"Ball Arrangement: {grid[0]} x {grid[1]} x {grid[2]}")
        print(f"Total Energy Balls: {n_balls}")
        print(f"Equation for Ball Cost: {n_balls} balls * {COST_PER_BALL_USD} usd/ball = {int(cost_balls)} usd")
        
        print("\nContainer Details:")
        print(f"Dimensions (cm): Length={dims[0]}, Width={dims[1]}, Height={dims[2]}")
        print(f"Equation for Surface Area: 2 * ({dims[0]}*{dims[1]} + {dims[1]}*{dims[2]} + {dims[2]}*{dims[0]}) = {int(surface_area)} cm^2")
        print(f"Equation for Material Cost: {int(surface_area)} cm^2 * {MATERIAL_COST_PER_CM2_USD} usd/cm^2 = {int(cost_material)} usd")

        print("\n--- Final Cost Calculation ---")
        print(f"Total Cost C = Ball Cost + Material Cost")
        print(f"C = {int(cost_balls)} + {int(cost_material)} = {int(total_cost)}")
        
        # Final answer format
        print(f"\n<<<{int(total_cost)}>>>")

    else:
        print("No solution found.")
        print("<<<0>>>")

solve_container_problem()