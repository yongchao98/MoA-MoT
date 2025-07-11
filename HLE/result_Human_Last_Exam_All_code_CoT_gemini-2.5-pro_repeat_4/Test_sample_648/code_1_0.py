import math

def solve_container_problem():
    """
    Solves the container design optimization problem.
    """
    # --- Problem Parameters ---
    ENERGY_TOTAL_MIN = 1000  # MJ
    ENERGY_PER_BALL = 30  # MJ
    COST_PER_BALL = 1000  # USD
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = BALL_RADIUS * 2.0  # cm
    MAX_SURFACE_AREA = 1000.0  # cm^2
    MATERIAL_COST_PER_CM2 = 200.0  # USD/cm^2
    PRECISION = 0.5  # cm

    # --- Step 1: Calculate the number of balls and their cost ---
    num_balls_needed = math.ceil(ENERGY_TOTAL_MIN / ENERGY_PER_BALL)
    cost_of_balls = num_balls_needed * COST_PER_BALL

    min_total_cost = float('inf')
    best_design = None

    # --- Step 2: Analyze Box Container ---
    # To hold N balls, we need nx*ny*nz >= N balls.
    # The minimal surface area for a given (nx, ny, nz) is when
    # L=nx*d, W=ny*d, H=nz*d.
    # SA_min = 2 * (d*d*nx*ny + d*d*ny*nz + d*d*nz*nx)
    # SA_min = 2 * d^2 * (nx*ny + ny*nz + nz*nx)
    # The most "cubical" arrangement of factors minimizes the surface area.
    # For 34 balls, the closest integer factors are (3, 3, 4) [product=36].
    nx, ny, nz = 3, 3, 4
    min_box_sa = 2 * (BALL_DIAMETER**2) * (nx*ny + ny*nz + nz*nx)
    # min_box_sa = 2 * 16 * (9 + 12 + 12) = 32 * 33 = 1056 cm^2
    # Since 1056 > 1000, no simple grid-packed box is feasible.
    # We will proceed assuming a box is not a valid solution.

    # --- Step 3: Analyze Cylinder Container ---
    # We check designs based on the number of layers of balls.
    # R/r ratios for packing k circles in a circle (from literature)
    # We need R_needed = r_ball * ratio = 2 * ratio
    k_min_radii_cm = {
        1: 2.0, 2: 4.0, 3: 4.309, 4: 4.828, 5: 5.464, 6: 6.0,
        7: 6.0, 8: 6.609, 9: 6.0, 10: 6.748, 11: 7.056, 12: 7.309,
        13: 7.513, 14: 7.786, 15: 8.053, 16: 8.163, 17: 8.426, 18: 8.532,
        19: 8.532, 20: 8.944, 34: 10.243
    }
    
    # We need a more precise value for k=17 (R/r = 3.813)
    k_min_radii_cm[17] = 2.0 * 3.813

    for num_layers in range(1, 10):
        # Calculate required height H
        min_h = num_layers * BALL_DIAMETER
        # Adjust H to precision
        h = math.ceil(min_h / PRECISION) * PRECISION

        # Calculate number of balls needed per layer
        k_per_layer = math.ceil(num_balls_needed / num_layers)
        if k_per_layer not in k_min_radii_cm:
            continue # Skip if we don't have packing data

        # Calculate required radius R
        min_r = k_min_radii_cm[k_per_layer]
        # Adjust R to precision
        r = math.ceil(min_r / PRECISION) * PRECISION

        # Calculate surface area and cost
        surface_area = 2 * math.pi * r**2 + 2 * math.pi * r * h
        if surface_area <= MAX_SURFACE_AREA:
            container_cost = surface_area * MATERIAL_COST_PER_CM2
            total_cost = container_cost + cost_of_balls
            
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design = {
                    "type": "Cylinder",
                    "R": r,
                    "H": h,
                    "SA": surface_area,
                    "cost": total_cost,
                    "balls": num_balls_needed,
                    "cost_balls": cost_of_balls,
                    "cost_container": container_cost
                }

    # --- Step 4: Output the result ---
    if best_design:
        C = round(best_design['cost'])
        SA_val = best_design['SA']
        cost_mat = MATERIAL_COST_PER_CM2
        n_balls = best_design['balls']
        cost_ball = COST_PER_BALL
        
        print("Optimal Design Found: Cylinder")
        print(f"Dimensions: Radius = {best_design['R']:.1f} cm, Height = {best_design['H']:.1f} cm")
        print(f"This design holds {best_design['balls']} balls in {int(best_design['H']/BALL_DIAMETER)} layer(s).")
        print("\nCost Calculation:")
        # The user wants each number in the final equation printed.
        print(f"Total Cost = (Surface Area * Material Cost) + (Number of Balls * Ball Cost)")
        print(f"C = ({SA_val:.2f} cm^2 * {cost_mat:.0f} USD/cm^2) + ({n_balls} * {cost_ball} USD)")
        
        container_cost_val = best_design['cost_container']
        ball_cost_val = best_design['cost_balls']
        print(f"C = {container_cost_val:.2f} USD + {ball_cost_val:.2f} USD")
        
        print(f"C = {C}")
        print(f"<<<{C}>>>")

    else:
        print("No feasible design found that meets all constraints.")
        print("<<<0>>>")

solve_container_problem()