import math

def solve_pioneer_probe_design():
    """
    Finds the minimum cost design for a container to hold energy balls.
    """
    # Problem Constants
    ENERGY_PER_BALL_MJ = 30
    REQUIRED_ENERGY_MJ = 1000
    BALL_COST_USD = 1000
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4

    MATERIAL_COST_PER_CM2 = 200
    MAX_SURFACE_AREA_CM2 = 1000
    PRECISION_CM = 0.5

    # Step 1: Calculate minimum number of balls required
    min_balls_needed = math.ceil(REQUIRED_ENERGY_MJ / ENERGY_PER_BALL_MJ)
    
    # Initialize best solution tracker
    best_design = {'cost': float('inf')}

    # --- Helper function for circle packing in a cylinder layer ---
    def get_min_cylinder_radius(balls_in_layer):
        # Pre-calculated minimum radii for packing n circles of radius r_ball=2
        # (Based on known circle packing results, r_packing = k * r_ball)
        r_ball = BALL_RADIUS_CM
        min_radii = {
            1: 1 * r_ball, 2: 2 * r_ball, 3: (1 + 2 / math.sqrt(3)) * r_ball,
            4: (1 + math.sqrt(2)) * r_ball, 5: 2.701 * r_ball, # Approximate k-value
            6: 2 * r_ball, # This is for a specific compact arrangement, not a hexagon row
            7: 3 * r_ball
        }
        if balls_in_layer > 7:
            # For larger n, approximate as a dense pack
            # Area_balls / packing_density = Area_circle
            # n * pi * r_b^2 / 0.9069 = pi * R^2 => R = r_b * sqrt(n/0.9069)
            return r_ball * math.sqrt(balls_in_layer / 0.9069)
        
        return min_radii.get(balls_in_layer, float('inf'))

    # --- Function to find the best cylinder design ---
    def find_best_cylinder():
        best_cylinder = {'cost': float('inf')}
        for n_layer in range(1, 20): # Iterate through possible number of balls per layer
            r_min_ideal = get_min_cylinder_radius(n_layer)
            # Enforce manufacturing precision
            R = math.ceil(r_min_ideal / PRECISION_CM) * PRECISION_CM
            
            n_z = math.ceil(min_balls_needed / n_layer)
            H = n_z * BALL_DIAMETER_CM

            SA = 2 * math.pi * R**2 + 2 * math.pi * R * H
            
            if SA <= MAX_SURFACE_AREA_CM2:
                cost_container = SA * MATERIAL_COST_PER_CM2
                cost_balls = min_balls_needed * BALL_COST_USD
                total_cost = cost_container + cost_balls
                
                if total_cost < best_cylinder['cost']:
                    best_cylinder = {
                        'type': 'Cylinder', 'cost': total_cost, 'R': R, 'H': H, 'SA': SA,
                        'layout': (n_layer, n_z), 'balls_fit': n_layer * n_z
                    }
        return best_cylinder if best_cylinder['cost'] != float('inf') else None

    # --- Function to find the best box design ---
    def find_best_box():
        best_box = {'cost': float('inf')}
        max_dim_steps = int(40 / PRECISION_CM) # Search dimensions up to 40cm
        
        for h_steps in range(1, max_dim_steps + 1):
            H = h_steps * PRECISION_CM
            for w_steps in range(h_steps, max_dim_steps + 1):
                W = w_steps * PRECISION_CM
                for l_steps in range(w_steps, max_dim_steps + 1):
                    L = l_steps * PRECISION_CM
                    
                    SA = 2 * (L*W + W*H + H*L)
                    if SA > MAX_SURFACE_AREA_CM2:
                        break # Increasing L further won't help

                    balls_fit_L = math.floor(L / BALL_DIAMETER_CM)
                    balls_fit_W = math.floor(W / BALL_DIAMETER_CM)
                    balls_fit_H = math.floor(H / BALL_DIAMETER_CM)
                    balls_fit = balls_fit_L * balls_fit_W * balls_fit_H

                    if balls_fit >= min_balls_needed:
                        cost_container = SA * MATERIAL_COST_PER_CM2
                        cost_balls = min_balls_needed * BALL_COST_USD
                        total_cost = cost_container + cost_balls
                        if total_cost < best_box['cost']:
                             best_box = {
                                'type': 'Box', 'cost': total_cost, 'L': L, 'W': W, 'H': H, 'SA': SA,
                                'balls_fit': balls_fit
                            }
        return best_box if best_box['cost'] != float('inf') else None


    # Step 2: Evaluate both container types
    best_cylinder_design = find_best_cylinder()
    best_box_design = find_best_box()

    # Step 3: Select the best overall design
    if best_cylinder_design and best_cylinder_design['cost'] < best_design['cost']:
        best_design = best_cylinder_design
    if best_box_design and best_box_design['cost'] < best_design['cost']:
        best_design = best_box_design
        
    # Step 4: Output the result
    if best_design['cost'] == float('inf'):
        print("No valid container design found that meets all criteria.")
        final_cost = 0
    else:
        print("Optimal Design Found:")
        print(f"Container Type: {best_design['type']}")
        if best_design['type'] == 'Cylinder':
             print(f"Dimensions: Radius = {best_design['R']:.1f} cm, Height = {best_design['H']:.1f} cm")
             print(f"Layout: {best_design['balls_fit']} balls ({best_design['layout'][0]} per layer, {best_design['layout'][1]} layers)")
        else: # Box
             print(f"Dimensions: L={best_design['L']:.1f}, W={best_design['W']:.1f}, H={best_design['H']:.1f} cm")
             print(f"Layout: {best_design['balls_fit']} balls")
        print("-" * 20)
        
        # Recalculate and print final equation for clarity
        sa = best_design['SA']
        cost_container = sa * MATERIAL_COST_PER_CM2
        cost_balls = min_balls_needed * BALL_COST_USD
        total_cost = cost_container + cost_balls
        final_cost = round(total_cost)

        print("Final Cost Calculation:")
        print(f"Minimum balls to buy for {REQUIRED_ENERGY_MJ} MJ = {min_balls_needed}")
        print(f"Cost of balls = {min_balls_needed} * ${BALL_COST_USD} = ${cost_balls}")
        print(f"Container surface area = {sa:.2f} cm2")
        print(f"Cost of container material = {sa:.2f} cm2 * ${MATERIAL_COST_PER_CM2}/cm2 = ${cost_container:.2f}")
        print(f"Total cost = ${cost_balls} + ${cost_container:.2f} = ${total_cost:.2f}")
        
    print(f"\nFinal Answer (total cost C): {final_cost}")
    print(f"<<<{final_cost}>>>")

solve_pioneer_probe_design()