import math

def solve():
    """
    Solves the container optimization problem.
    """
    # Problem constants
    ENERGY_PER_BALL_MJ = 25
    REQUIRED_ENERGY_MJ = 1000
    COST_PER_BALL_USD = 1000
    COST_PER_CM2_USD = 200
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    PRECISION_CM = 0.5

    min_balls = math.ceil(REQUIRED_ENERGY_MJ / ENERGY_PER_BALL_MJ)
    
    min_total_cost = float('inf')
    best_design_details = {}

    # --- Step 1: Analyze Box Container ---
    # For a box, using more than the minimum number of balls (40) will always increase the cost.
    # We find the best arrangement for N=40.
    N_box = min_balls
    best_box_area = float('inf')
    best_box_dims = (0, 0, 0)

    for nx in range(1, N_box + 1):
        if N_box % nx == 0:
            for ny in range(1, (N_box // nx) + 1):
                if (N_box // nx) % ny == 0:
                    nz = (N_box // nx) // ny
                    l = nx * BALL_DIAMETER_CM
                    w = ny * BALL_DIAMETER_CM
                    h = nz * BALL_DIAMETER_CM
                    area = 2 * (l*w + l*h + w*h)
                    if area < best_box_area:
                        best_box_area = area
                        best_box_dims = (l, w, h)

    cost_balls_box = N_box * COST_PER_BALL_USD
    cost_container_box = best_box_area * COST_PER_CM2_USD
    total_cost_box = cost_balls_box + cost_container_box

    if total_cost_box < min_total_cost:
        min_total_cost = total_cost_box
        best_design_details = {
            "type": "Box", "N": N_box, "area": best_box_area,
            "cost": total_cost_box, "dims": best_box_dims
        }

    # --- Step 2: Analyze Cylinder Container ---
    # For a cylinder, a more efficient packing might be achieved with N > 40.
    # We search for N in a reasonable range (40 to 100).

    # Pre-computed optimal packing ratios (R_container / r_circle) for n circles in a circle.
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cinc/cinc.html
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.154, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0, 8: 3.304, 
        9: 3.512, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236, 14: 4.328, 
        15: 4.521, 16: 4.615, 17: 4.792, 18: 4.863, 19: 5.0, 20: 5.122,
        21: 5.248, 22: 5.343, 23: 5.517, 24: 5.607, 25: 5.753, 26: 5.861,
        27: 5.968, 28: 6.079, 29: 6.136, 30: 6.257
    }

    for n_balls in range(min_balls, 101):
        for n_layers in range(1, n_balls + 1):
            if n_balls % n_layers == 0:
                balls_per_layer = n_balls // n_layers
                if balls_per_layer in packing_ratios:
                    # Calculate cylinder dimensions
                    height = n_layers * BALL_DIAMETER_CM
                    
                    # Radius must contain the packed circles and meet precision
                    min_radius_req = packing_ratios[balls_per_layer] * BALL_RADIUS_CM
                    radius = math.ceil(min_radius_req / PRECISION_CM) * PRECISION_CM

                    # Calculate area and cost
                    area = (2 * math.pi * radius**2) + (2 * math.pi * radius * height)
                    cost_balls = n_balls * COST_PER_BALL_USD
                    cost_container = area * COST_PER_CM2_USD
                    total_cost = cost_balls + cost_container

                    if total_cost < min_total_cost:
                        min_total_cost = total_cost
                        best_design_details = {
                            "type": "Cylinder", "N": n_balls, "area": area,
                            "cost": total_cost, "dims": (radius, height),
                            "layout": (n_layers, balls_per_layer)
                        }

    # --- Step 3: Print the final answer ---
    if best_design_details:
        design = best_design_details["type"]
        N = best_design_details["N"]
        area = best_design_details["area"]
        cost = best_design_details["cost"]
        ball_cost = N * COST_PER_BALL_USD
        container_cost = area * COST_PER_CM2_USD
        
        print(f"The optimal design is a {design} container.")
        print(f"It holds {N} energy balls.")

        if design == "Box":
            l, w, h = best_design_details["dims"]
            print(f"Dimensions (L, W, H): {l} cm, {w} cm, {h} cm.")
            print("\n--- Total Cost Calculation ---")
            print(f"Total Cost = (N_balls * Cost_per_ball) + (Surface_Area * Cost_per_cm2)")
            print(f"Total Cost = ({N} * {COST_PER_BALL_USD}) + (2 * ({l}*{w} + {l}*{h} + {w}*{h}) * {COST_PER_CM2_USD})")
            print(f"Total Cost = {ball_cost} + ({area:.2f} * {COST_PER_CM2_USD})")
            print(f"Total Cost = {ball_cost} + {container_cost:.2f}")

        elif design == "Cylinder":
            R, H = best_design_details["dims"]
            print(f"Cylinder Radius: {R} cm, Cylinder Height: {H} cm.")
            print(f"Arrangement: {best_design_details['layout'][0]} layers of {best_design_details['layout'][1]} balls.")
            print("\n--- Total Cost Calculation ---")
            print(f"Total Cost = (N_balls * Cost_per_ball) + (Surface_Area * Cost_per_cm2)")
            print(f"Total Cost = ({N} * {COST_PER_BALL_USD}) + (2 * pi * R^2 + 2 * pi * R * H) * {COST_PER_CM2_USD}")
            print(f"Total Cost = ({N} * {COST_PER_BALL_USD}) + (2 * {math.pi:.4f} * {R}^2 + 2 * {math.pi:.4f} * {R} * {H}) * {COST_PER_CM2_USD}")
            print(f"Total Cost = {ball_cost} + ({area:.2f} * {COST_PER_CM2_USD})")
            print(f"Total Cost = {ball_cost} + {container_cost:.2f}")
        
        print(f"\nFinal Total Cost = {cost:.2f} USD")
        print(f"\n<<<{cost:.2f}>>>")

    else:
        print("No solution found.")
        print("\n<<<0>>>")

if __name__ == '__main__':
    solve()