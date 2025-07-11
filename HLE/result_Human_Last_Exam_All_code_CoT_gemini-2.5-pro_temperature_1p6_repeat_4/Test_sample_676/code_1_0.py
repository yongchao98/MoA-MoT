import math

def solve_container_problem():
    """
    Solves the container design optimization problem.
    """

    # --- Problem Constants ---
    ENERGY_PER_BALL = 25  # MJ
    COST_PER_BALL = 1000  # USD
    MATERIAL_COST_PER_CM2 = 200  # USD
    TOTAL_ENERGY_REQUIRED = 1000  # MJ
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = 4.0  # cm
    PRECISION = 0.5  # cm

    # --- Step 1: Calculate minimum number of energy balls ---
    N_BALLS_MIN = math.ceil(TOTAL_ENERGY_REQUIRED / ENERGY_PER_BALL)

    # --- Step 2: Box Container Optimization ---
    min_cost_box = float('inf')
    best_box_config = {}

    # Iterate through possible numbers of balls along each dimension (nx, ny, nz)
    limit = N_BALLS_MIN + 5 # Search up to a reasonable limit
    for nx in range(1, limit):
        for ny in range(nx, limit):
            if nx * ny > N_BALLS_MIN * 2 and ny > nx * 2: # Pruning the search space
                break
            
            nz = math.ceil(N_BALLS_MIN / (nx * ny))
            n_actual = nx * ny * nz

            # Container dimensions are determined by the number of balls.
            L = nx * BALL_DIAMETER
            W = ny * BALL_DIAMETER
            H = nz * BALL_DIAMETER
            
            # Calculate surface area and total cost
            surface_area = 2 * (L * W + L * H + W * H)
            cost_balls = n_actual * COST_PER_BALL
            cost_material = surface_area * MATERIAL_COST_PER_CM2
            total_cost = cost_balls + cost_material
            
            if total_cost < min_cost_box:
                min_cost_box = total_cost
                best_box_config = {
                    'type': 'Box', 'cost': total_cost, 'balls_arrangement': (nx, ny, nz),
                    'balls_total': n_actual, 'dimensions_cm': (L, W, H),
                    'surface_area_cm2': surface_area, 'cost_balls': cost_balls, 'cost_material': cost_material
                }

    # --- Step 3: Cylinder Container Optimization ---
    min_cost_cyl = float('inf')
    best_cyl_config = {}

    # Ratio of container radius to ball radius for packing k balls in a layer.
    R_k_ratios = {
        1: 1.0, 2: 2.0, 3: 2.1547, 4: 2.4142, 5: 2.7013, 6: 3.0, 7: 3.0,
        8: 3.3048, 9: 3.5129, 10: 3.8130, 11: 3.9238, 12: 4.0296, 13: 4.2361,
        14: 4.3284, 15: 4.5213, 16: 4.6155, 17: 4.7923, 18: 4.8637, 19: 4.8637,
        20: 5.1223, 21: 5.2058, 22: 5.3229, 23: 5.5342, 24: 5.6063, 25: 5.6883,
        26: 5.8770, 27: 6.0, 28: 6.0594, 29: 6.1547, 30: 6.2730, 31: 6.3400,
        32: 6.4678, 33: 6.5383, 34: 6.6234, 35: 6.7313, 36: 6.8394, 37: 6.8672,
        38: 7.0270, 39: 7.0860, 40: 7.1754,
    }

    # Iterate through possible numbers of balls per layer (k)
    for k in range(1, N_BALLS_MIN + 1):
        m = math.ceil(N_BALLS_MIN / k) # Number of layers
        n_actual = m * k

        H = m * BALL_DIAMETER
        # Radius must be rounded UP to the next multiple of the precision
        ideal_R = R_k_ratios[k] * BALL_RADIUS
        R = math.ceil(ideal_R / PRECISION) * PRECISION

        surface_area = (2 * math.pi * R**2) + (2 * math.pi * R * H)
        cost_balls = n_actual * COST_PER_BALL
        cost_material = surface_area * MATERIAL_COST_PER_CM2
        total_cost = cost_balls + cost_material

        if total_cost < min_cost_cyl:
            min_cost_cyl = total_cost
            best_cyl_config = {
                'type': 'Cylinder', 'cost': total_cost, 'balls_arrangement': (k, m),
                'balls_total': n_actual, 'dimensions_cm': (R, H),
                'surface_area_cm2': surface_area, 'cost_balls': cost_balls, 'cost_material': cost_material
            }

    # --- Step 4: Compare and Final Answer ---
    if min_cost_box < min_cost_cyl:
        final_design = best_box_config
    else:
        final_design = best_cyl_config

    C = final_design['cost']

    print("Pioneer Probe Container Design Optimization\n")
    print(f"1. Energy Requirement: >= {TOTAL_ENERGY_REQUIRED} MJ")
    print(f"   - Minimum balls needed: {N_BALLS_MIN}\n")

    print(f"2. Optimal Design Found: {final_design['type']}\n")
    
    print("3. Final Design Details:")
    if final_design['type'] == 'Box':
        L, W, H = final_design['dimensions_cm']
        print(f"   - Ball Arrangement (in layers): {final_design['balls_arrangement'][0]}x{final_design['balls_arrangement'][1]}x{final_design['balls_arrangement'][2]}")
        print(f"   - Container Dimensions (L x W x H): {L}cm x {W}cm x {H}cm")
    else: # Cylinder
        R, H = final_design['dimensions_cm']
        print(f"   - Ball Arrangement: {final_design['balls_arrangement'][1]} layers of {final_design['balls_arrangement'][0]} balls")
        print(f"   - Container Dimensions (Radius x Height): {R}cm x {H}cm")
    print(f"   - Total balls used: {final_design['balls_total']}")
    print(f"   - Container surface area: {final_design['surface_area_cm2']:.2f} cm^2\n")

    print("4. Cost Calculation:")
    cost_of_balls = final_design['balls_total'] * COST_PER_BALL
    cost_of_material = final_design['cost'] - cost_of_balls
    
    print("   Total Cost C = (Total Balls * Cost per Ball) + (Surface Area * Material Cost per cm^2)\n")
    
    # Outputting each number in the final equation
    if final_design['type'] == 'Box':
        L, W, H = final_design['dimensions_cm']
        print(f"   C = ({final_design['balls_total']} * {COST_PER_BALL}) + (2 * ({L}*{W} + {L}*{H} + {W}*{H}) * {MATERIAL_COST_PER_CM2})")
        print(f"   C = {final_design['cost_balls']:,.2f} + ({final_design['surface_area_cm2']:.4f} * {MATERIAL_COST_PER_CM2})")
    else: # Cylinder
        R, H = final_design['dimensions_cm']
        print(f"   C = ({final_design['balls_total']} * {COST_PER_BALL}) + ((2 * pi * {R}^2 + 2 * pi * {R} * {H}) * {MATERIAL_COST_PER_CM2})")
        print(f"   C = {final_design['cost_balls']:,.2f} + ({final_design['surface_area_cm2']:.4f} * {MATERIAL_COST_PER_CM2})")

    print(f"   C = ${cost_of_balls:,.2f} + ${cost_of_material:,.2f}")
    print(f"   Total Cost = ${C:,.2f}")
    
    print(f"\n<<<{int(round(C))}>>>")

if __name__ == '__main__':
    solve_container_problem()