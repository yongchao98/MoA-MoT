import math

def solve_pioneer_probe_design():
    """
    Solves the container design optimization problem for the Pioneer probe.
    It analyzes both box and cylinder containers and finds the one with the
    lowest total cost for materials and energy balls.
    """
    # --- Problem Parameters ---
    ENERGY_REQ_MJ = 1000
    ENERGY_PER_BALL_MJ = 30
    BALL_COST_USD = 1000
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    MAX_SURFACE_AREA_CM2 = 1000
    MATERIAL_COST_PER_CM2 = 200
    PRECISION_CM = 0.5

    # --- Step 1: Minimum number of balls ---
    min_balls_needed = math.ceil(ENERGY_REQ_MJ / ENERGY_PER_BALL_MJ)

    min_total_cost = float('inf')
    best_design = None

    # --- Step 2: Analyze Box Container ---
    # We will search for a valid box design by iterating through possible
    # numbers of balls arranged in a grid (n_L x n_W x n_H).
    # Based on mathematical analysis (AM-GM inequality on surface area vs volume),
    # the minimum surface area for a box holding >= 34 balls is > 1000 cm^2.
    # For example, for a 4x3x3 arrangement (36 balls), the tightest box is
    # 16x12x12 cm. Area = 2*(16*12 + 16*12 + 12*12) = 1056 cm^2.
    # No denser packing of >= 34 balls results in a smaller surface area.
    # Therefore, the box container is not a feasible option.

    # --- Step 3 & 4: Analyze Cylinder Container & Search for Optimal Solution ---
    # Data for packing n circles in a larger circle.
    # k_values[n] = R/r, where R is the radius of the container circle and
    # r is the radius of the packed circles.
    # Source: packomania.com
    k_values = {
        1: 1.0, 2: 2.0, 3: 2.1548, 4: 2.4143, 5: 2.7013, 6: 3.0, 7: 3.0,
        8: 3.3048, 9: 3.6131, 10: 3.8130, 11: 3.9715, 12: 4.0296, 13: 4.2361,
        14: 4.3284, 15: 4.5210, 16: 4.6154, 17: 4.7924, 18: 4.8637, 19: 4.8637,
        20: 5.1223
    }
    
    # Iterate through number of layers (n_H) and balls per layer (n_layer)
    for n_H in range(1, min_balls_needed + 1):
        for n_layer, k in k_values.items():
            num_balls = n_H * n_layer
            
            if num_balls < min_balls_needed:
                continue

            # Calculate container dimensions
            height = n_H * BALL_DIAMETER_CM
            
            # Required radius for packing, then adjusted for precision
            min_radius_pack = k * BALL_RADIUS_CM
            container_radius = math.ceil(min_radius_pack / PRECISION_CM) * PRECISION_CM

            # Calculate surface area
            surface_area = 2 * math.pi * container_radius * (container_radius + height)
            
            if surface_area <= MAX_SURFACE_AREA_CM2:
                # Calculate total cost
                cost_balls = num_balls * BALL_COST_USD
                cost_container = surface_area * MATERIAL_COST_PER_CM2
                total_cost = cost_balls + cost_container

                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        "type": "Cylinder",
                        "num_balls": num_balls,
                        "n_H": n_H,
                        "n_layer": n_layer,
                        "height": height,
                        "radius": container_radius,
                        "surface_area": surface_area,
                        "cost": total_cost,
                        "cost_balls": cost_balls,
                        "cost_container": cost_container
                    }

    # --- Step 5: Output the result ---
    if best_design:
        print("Optimal design found:")
        print(f"Container Type: {best_design['type']}")
        print(f"Configuration: {best_design['n_H']} layers of {best_design['n_layer']} balls")
        print(f"Dimensions (H x R): {best_design['height']:.1f} cm x {best_design['radius']:.1f} cm")
        print("\nFinal Cost Calculation:")
        print(f"C = (Number of Balls * Cost per Ball) + (Surface Area * Cost per cm^2)")
        
        N_val = best_design['num_balls']
        N_cost = BALL_COST_USD
        A_val = best_design['surface_area']
        A_cost = MATERIAL_COST_PER_CM2
        C_val = best_design['cost']
        
        print(f"C = {N_val} * {N_cost} + {A_val:.2f} * {A_cost} = {C_val:.2f}")

        # The final answer in the required format
        print(f"\n<<<{C_val:.2f}>>>")

    else:
        print("No feasible design solution found that meets all constraints.")
        print("<<<0>>>")

solve_pioneer_probe_design()