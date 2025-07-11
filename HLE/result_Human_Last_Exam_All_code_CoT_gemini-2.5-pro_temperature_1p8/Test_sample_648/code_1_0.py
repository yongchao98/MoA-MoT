import math
import sys

def solve_pioneer_probe_design():
    """
    Solves the Pioneer probe container design problem by systematically searching
    for the lowest-cost configuration that meets all constraints.
    """
    
    # --- Problem Constants ---
    BALL_DIAMETER = 4.0  # cm (radius is 2 cm)
    BALL_ENERGY = 30     # MJ
    BALL_COST = 1000     # USD
    
    MATERIAL_COST_PER_CM2 = 200  # USD/cm^2
    MAX_SURFACE_AREA = 1000      # cm^2
    
    REQUIRED_ENERGY = 1000       # MJ
    PRECISION = 0.5              # cm

    # --- Step 1: Calculate minimum number of balls ---
    # The total cost function increases with the number of balls, so we only
    # need to find a solution for the minimum required number of balls.
    n_balls = math.ceil(REQUIRED_ENERGY / BALL_ENERGY)
    
    min_total_cost = float('inf')
    best_design = None
    cost_balls = n_balls * BALL_COST

    # --- Step 2: Search for the best Box Container ---
    # To minimize surface area for a given volume, dimensions should be as
    # close as possible (cube-like). We iterate nx, ny, nz such that their
    # product is >= n_balls, and check the resulting surface area.
    for nx in range(1, n_balls + 1):
        for ny in range(nx, n_balls + 1):
            if nx * ny > n_balls:
                break
            nz = math.ceil(n_balls / (nx * ny))
            
            L = nx * BALL_DIAMETER
            W = ny * BALL_DIAMETER
            H = nz * BALL_DIAMETER
            
            sa = 2 * (L * W + W * H + H * L)
            
            if sa <= MAX_SURFACE_AREA:
                total_cost = sa * MATERIAL_COST_PER_CM2 + cost_balls
                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        "type": "box", "cost": total_cost, "sa": sa, 
                        "n_balls": n_balls, "L": L, "W": W, "H": H
                    }
    
    # --- Step 3: Search for the best Cylinder Container ---
    # Iterate through possible numbers of layers for the balls.
    for n_layers in range(1, n_balls + 1):
        H = n_layers * BALL_DIAMETER
        
        # Balls needed per layer.
        n_per_layer = math.ceil(n_balls / n_layers)
        
        # Find the minimum radius to pack n_per_layer balls in a grid.
        min_R_for_layer = float('inf')
        best_nx_l = 0
        for nx_l in range(1, n_per_layer + 1):
            ny_l = math.ceil(n_per_layer / nx_l)
            
            # The required diameter is the diagonal of the bounding box (nx*d x ny*d).
            required_diameter = math.sqrt((nx_l * BALL_DIAMETER)**2 + (ny_l * BALL_DIAMETER)**2)
            required_radius = required_diameter / 2.0
            
            # The actual radius must be a multiple of the manufacturing precision.
            actual_radius = math.ceil(required_radius / PRECISION) * PRECISION
            
            if actual_radius < min_R_for_layer:
                min_R_for_layer = actual_radius
                
        R = min_R_for_layer
        sa = 2 * math.pi * R**2 + 2 * math.pi * R * H
        
        if sa <= MAX_SURFACE_AREA:
            total_cost = sa * MATERIAL_COST_PER_CM2 + cost_balls
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design = {
                    "type": "cylinder", "cost": total_cost, "sa": sa,
                    "n_balls": n_balls, "R": R, "H": H
                }

    # --- Step 4: Output the result ---
    if best_design:
        cost_of_material = best_design['sa'] * MATERIAL_COST_PER_CM2
        cost_of_balls = best_design['n_balls'] * BALL_COST
        
        print("A feasible design was found.")
        print(f"Optimal design is a {best_design['type']}.")
        
        print("\n--- Cost Calculation ---")
        print("Total Cost = (Surface Area * Material Cost/cm^2) + (Number of Balls * Ball Cost)")
        # Outputting each number in the final equation as requested
        print(f"C = {best_design['sa']:.2f} * {MATERIAL_COST_PER_CM2} + {best_design['n_balls']} * {BALL_COST}")
        print(f"C = {cost_of_material:.2f} + {cost_of_balls}")
        final_answer = f"{best_design['cost']:.0f}"

    else:
        print("No feasible design was found that meets all constraints (especially the max surface area of 1000 cm^2).")
        final_answer = "0"
        
    print(f"\n<<<{final_answer}>>>")

solve_pioneer_probe_design()