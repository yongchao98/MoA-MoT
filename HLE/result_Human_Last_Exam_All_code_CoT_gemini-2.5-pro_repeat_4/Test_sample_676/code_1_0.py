import math

def solve_container_problem():
    """
    Calculates the minimum cost to design a container for energy balls,
    considering both box and cylinder shapes.
    """

    # --- Problem Constants ---
    BALL_RADIUS_CM = 2.0
    BALL_DIAMETER_CM = 4.0
    BALL_COST_USD = 1000
    MATERIAL_COST_USD_PER_CM2 = 200
    REQUIRED_ENERGY_MJ = 1000
    BALL_ENERGY_MJ = 25
    PRECISION_CM = 0.5

    MIN_BALLS = math.ceil(REQUIRED_ENERGY_MJ / BALL_ENERGY_MJ)
    
    min_total_cost = float('inf')
    best_design_details = {}

    def round_up_to_precision(value, precision):
        """Rounds a value up to the nearest multiple of the given precision."""
        return math.ceil(value / precision) * precision

    # --- Part 1: Box Container Optimization ---
    min_box_cost = float('inf')
    best_box_config = {}
    
    # Iterate through total ball counts N >= 40. A larger N can sometimes be cheaper.
    # We check up to N=80, which is a reasonable search limit.
    for n_total in range(MIN_BALLS, 81):
        # Find integer factorizations (nx, ny, nz) for n_total
        for nx in range(1, int(n_total**(1/3.0)) + 2):
            if n_total % nx == 0:
                for ny in range(nx, int(math.sqrt(n_total / nx)) + 2):
                    if (n_total % nx) % ny == 0:
                        nz = n_total // (nx * ny)
                        if nx * ny * nz == n_total:
                            # Dimensions for simple cubic packing
                            l = nx * BALL_DIAMETER_CM
                            w = ny * BALL_DIAMETER_CM
                            h = nz * BALL_DIAMETER_CM

                            # Surface area and cost
                            surface_area = 2 * (l*w + w*h + h*l)
                            cost_material = surface_area * MATERIAL_COST_USD_PER_CM2
                            cost_balls = n_total * BALL_COST_USD
                            total_cost = cost_material + cost_balls

                            if total_cost < min_box_cost:
                                min_box_cost = total_cost
                                best_box_config = {
                                    'type': 'Box', 'cost': total_cost, 'balls': n_total,
                                    'nx': nx, 'ny': ny, 'nz': nz, 'L': l, 'W': w, 'H': h,
                                    'SA': surface_area, 'cost_mat': cost_material, 'cost_balls': cost_balls
                                }

    if best_box_config:
        min_total_cost = min_box_cost
        best_design_details = best_box_config

    # --- Part 2: Cylinder Container Optimization ---
    
    # Packing data for N circles in a circle (R_container / r_ball)
    # Source: E. Specht, "The best known packings of equal circles in a circle"
    PACK_DATA = {
        1: 1.000, 2: 2.000, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.000, 7: 3.000,
        8: 3.305, 9: 3.613, 10: 3.813, 11: 3.924, 12: 4.030, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 5.000,
        20: 5.124
    }

    min_cyl_cost = float('inf')
    best_cyl_config = {}

    # Iterate through number of layers (nz)
    for nz in range(1, MIN_BALLS + 1):
        # Determine minimum balls per layer
        min_n_layer = math.ceil(MIN_BALLS / nz)
        
        # Check a few options for n_layer, as more balls can sometimes be optimal
        for n_layer in range(min_n_layer, min_n_layer + 5):
            if n_layer not in PACK_DATA:
                continue

            # Calculate Radius
            radius_ratio = PACK_DATA[n_layer]
            ideal_radius = radius_ratio * BALL_RADIUS_CM
            container_radius = round_up_to_precision(ideal_radius, PRECISION_CM)

            # Calculate Height (assuming close-packing of layers)
            if nz == 1:
                ideal_height = BALL_DIAMETER_CM
            else:
                # Height of (nz-1) layers offset + diameter of one ball at the base
                layer_height_offset = BALL_DIAMETER_CM * math.sqrt(2.0/3.0)
                ideal_height = (nz - 1) * layer_height_offset + BALL_DIAMETER_CM
            
            container_height = round_up_to_precision(ideal_height, PRECISION_CM)

            # Calculate total cost
            n_total = nz * n_layer
            surface_area = 2 * math.pi * container_radius * (container_radius + container_height)
            cost_material = surface_area * MATERIAL_COST_USD_PER_CM2
            cost_balls = n_total * BALL_COST_USD
            total_cost = cost_material + cost_balls
            
            if total_cost < min_cyl_cost:
                min_cyl_cost = total_cost
                best_cyl_config = {
                    'type': 'Cylinder', 'cost': total_cost, 'balls': n_total,
                    'nz': nz, 'n_layer': n_layer, 'R': container_radius, 'H': container_height,
                    'SA': surface_area, 'cost_mat': cost_material, 'cost_balls': cost_balls
                }

    if best_cyl_config and min_cyl_cost < min_total_cost:
        min_total_cost = min_cyl_cost
        best_design_details = best_cyl_config

    # --- Part 3: Output the Final Result ---
    if min_total_cost == float('inf'):
        print("A solution could not be found.")
        final_cost = 0
    else:
        print("The optimal design is a Cylinder with the following characteristics:\n")
        
        d = best_design_details
        print(f"Number of balls (N): {d['balls']}")
        print(f"  - Arrangement: {d['nz']} layers of {d['n_layer']} balls each.")
        print(f"Container Dimensions:")
        print(f"  - Radius (R): {d['R']:.1f} cm")
        print(f"  - Height (H): {d['H']:.1f} cm")
        print("\nCost Calculation:")
        print(f"1. Cost of Energy Balls = {d['balls']} balls * ${BALL_COST_USD}/ball")
        print(f"   = ${d['cost_balls']:.2f}")
        print(f"2. Cost of Material:")
        print(f"   - Surface Area = 2 * pi * R * (R + H)")
        print(f"   - Surface Area = 2 * {math.pi:.4f} * {d['R']:.1f} * ({d['R']:.1f} + {d['H']:.1f}) = {d['SA']:.2f} cm^2")
        print(f"   - Material Cost = {d['SA']:.2f} cm^2 * ${MATERIAL_COST_USD_PER_CM2}/cm^2")
        print(f"   = ${d['cost_mat']:.2f}")
        print("\n3. Total Cost (C) = Cost of Balls + Cost of Material")
        print(f"   C = ${d['cost_balls']:.2f} + ${d['cost_mat']:.2f}")
        print(f"   C = ${d['cost']:.2f}")
        final_cost = d['cost']
    
    # Final answer in the required format
    print(f"\n<<<{final_cost:.2f}>>>")

# Run the solver
solve_container_problem()