import math

def solve_container_problem():
    """
    This function solves the container optimization problem by finding the lowest
    cost design for both a box and a cylinder, and then selecting the best one.
    """
    
    # --- Problem Constants ---
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = 4.0 # cm
    BALL_ENERGY = 25    # MJ
    BALL_COST = 1000    # USD
    MATERIAL_COST_PER_CM2 = 200 # USD
    ENERGY_REQUIREMENT = 1000 # MJ
    PRECISION = 0.5 # cm
    
    MIN_BALLS_REQUIRED = math.ceil(ENERGY_REQUIREMENT / BALL_ENERGY)

    min_total_cost = float('inf')
    best_design = {}

    # --- Part 1: Box Optimization ---
    min_box_cost = float('inf')
    best_box_design = {}
    
    # Search space for ball counts along each dimension. A large dimension is inefficient.
    # Searching up to 41 is sufficient as 1x1x40 is a valid, though inefficient, baseline.
    for nx in range(1, MIN_BALLS_REQUIRED + 1):
        for ny in range(nx, MIN_BALLS_REQUIRED + 1):
            for nz in range(ny, MIN_BALLS_REQUIRED + 1):
                num_balls = nx * ny * nz
                
                if num_balls < MIN_BALLS_REQUIRED:
                    continue

                L, W, H = BALL_DIAMETER * nx, BALL_DIAMETER * ny, BALL_DIAMETER * nz
                
                surface_area = 2 * (L*W + L*H + W*H)
                material_cost = surface_area * MATERIAL_COST_PER_CM2
                ball_cost = num_balls * BALL_COST
                total_cost = material_cost + ball_cost
                
                if total_cost < min_box_cost:
                    min_box_cost = total_cost
                    best_box_design = {
                        "type": "Box",
                        "nx": nx, "ny": ny, "nz": nz,
                        "L": L, "W": W, "H": H,
                        "num_balls": num_balls,
                        "surface_area": surface_area,
                        "ball_cost": ball_cost,
                        "material_cost": material_cost,
                        "total_cost": total_cost,
                    }
    
    if min_box_cost < min_total_cost:
        min_total_cost = min_box_cost
        best_design = best_box_design

    # --- Part 2: Cylinder Optimization ---
    # We use optimal 2D circle packing data for the layers.
    # k[n] is the ratio R_container/r_ball for packing n circles. r_ball = 2.
    k_vals = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 5.0,
        20: 5.122, 21: 5.211, 22: 5.322, 23: 5.431, 24: 5.545, 25: 5.655,
        26: 5.753, 27: 5.867, 28: 5.961, 29: 6.043, 30: 6.155, 31: 6.233,
        32: 6.314, 33: 6.398, 34: 6.480, 35: 6.574, 36: 6.643, 37: 6.747,
        38: 6.814, 39: 6.892, 40: 6.994
    }

    min_cyl_cost = float('inf')
    best_cyl_design = {}

    for n_layer in range(1, MIN_BALLS_REQUIRED + 1):
        num_layers = math.ceil(MIN_BALLS_REQUIRED / n_layer)
        num_balls = n_layer * num_layers
        
        # Calculate container dimensions
        H = BALL_DIAMETER * num_layers
        
        # Radius of container must enclose n_layer balls
        k = k_vals.get(n_layer)
        if k is None: continue # Skip if no data for this n_layer
            
        r_needed = BALL_RADIUS * k
        # Enforce precision by rounding up to nearest multiple of 0.5
        R = math.ceil(r_needed / PRECISION) * PRECISION
        
        surface_area = 2 * math.pi * R**2 + 2 * math.pi * R * H
        material_cost = surface_area * MATERIAL_COST_PER_CM2
        ball_cost = num_balls * BALL_COST
        total_cost = material_cost + ball_cost
        
        if total_cost < min_cyl_cost:
            min_cyl_cost = total_cost
            best_cyl_design = {
                "type": "Cylinder",
                "n_layer": n_layer, "num_layers": num_layers,
                "R": R, "H": H,
                "num_balls": num_balls,
                "surface_area": surface_area,
                "ball_cost": ball_cost,
                "material_cost": material_cost,
                "total_cost": total_cost,
            }

    if min_cyl_cost < min_total_cost:
        min_total_cost = min_cyl_cost
        best_design = best_cyl_design

    # --- Part 3: Final Answer ---
    print("--- Optimal Design Found ---")
    if not best_design:
        print("No solution found.")
        final_cost = 0
    else:
        print(f"Container Type: {best_design['type']}")
        if best_design['type'] == 'Box':
            print(f"Dimensions (cm): L={best_design['L']}, W={best_design['W']}, H={best_design['H']}")
            print(f"Ball arrangement: {best_design['nx']}x{best_design['ny']}x{best_design['nz']}")
        else: # Cylinder
            print(f"Dimensions (cm): R={best_design['R']}, H={best_design['H']}")
            print(f"Ball arrangement: {best_design['num_layers']} layers of {best_design['n_layer']} balls")
        
        print(f"Total Balls: {best_design['num_balls']}")
        print(f"Container Surface Area: {best_design['surface_area']:.2f} cm^2")
        print("\n--- Cost Breakdown ---")
        
        # Outputting the numbers for the final equation as requested
        C = best_design['total_cost']
        cost_b = best_design['ball_cost']
        cost_m = best_design['material_cost']
        n_b = best_design['num_balls']
        p_b = BALL_COST
        area = best_design['surface_area']
        p_m = MATERIAL_COST_PER_CM2
        
        print(f"Final Equation: Total Cost = (Num Balls * Cost per Ball) + (Surface Area * Cost per cm^2)")
        print(f"Calculation: {C:.2f} = ({n_b} * {p_b}) + ({area:.2f} * {p_m})")
        print(f"Expanded:    {C:.2f} = {cost_b:.2f} + {cost_m:.2f}")

        final_cost = C

    return final_cost


# --- Execute the function and print the final answer in the required format ---
final_total_cost = solve_container_problem()
print(f"\n<<<C = {final_total_cost:.2f}>>>")
print(f"<<<{final_total_cost}>>>")
