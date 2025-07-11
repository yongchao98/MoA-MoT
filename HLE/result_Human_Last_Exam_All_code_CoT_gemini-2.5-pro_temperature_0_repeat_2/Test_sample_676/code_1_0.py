import math

def solve_container_problem():
    """
    Calculates the minimum cost to build a container for energy balls,
    considering both box and cylinder shapes.
    """

    # --- Problem Constants ---
    BALL_RADIUS = 2.0
    BALL_DIAMETER = 4.0
    BALL_ENERGY = 25  # MJ
    BALL_COST = 1000  # USD
    MATERIAL_COST_PER_CM2 = 200  # USD
    TOTAL_ENERGY_REQ = 1000  # MJ
    PRECISION = 0.5  # cm

    # --- Step 1: Calculate Ball Requirements ---
    NUM_BALLS = math.ceil(TOTAL_ENERGY_REQ / BALL_ENERGY)
    TOTAL_BALL_COST = NUM_BALLS * BALL_COST

    # --- Step 2: Optimize Box Container ---
    min_box_area = float('inf')
    best_box_config = None
    
    # We only need to check for N=40 balls, as packing more is less efficient.
    n_total = int(NUM_BALLS)
    # Iterate through all factor triplets of n_total
    for nx in range(1, n_total + 1):
        if n_total % nx == 0:
            for ny in range(nx, n_total // nx + 1):
                if (n_total % (nx * ny)) == 0:
                    nz = n_total // (nx * ny)
                    if nz >= ny:
                        # Dimensions based on simple grid packing
                        L = nx * BALL_DIAMETER
                        W = ny * BALL_DIAMETER
                        H = nz * BALL_DIAMETER
                        
                        area = 2 * (L*W + W*H + H*L)
                        if area < min_box_area:
                            min_box_area = area
                            best_box_config = (nx, ny, nz)

    box_material_cost = min_box_area * MATERIAL_COST_PER_CM2
    total_box_cost = box_material_cost + TOTAL_BALL_COST

    # --- Step 3: Optimize Cylinder Container ---
    # R/r ratios for packing N circles in a circle. Source: hydra.nat.uni-magdeburg.de/packing/cinc/
    packing_ratios = {
        1: 1.0, 2: 1.0, 3: 1.155, 4: 1.414, 5: 1.618, 6: 1.732, 7: 2.0,
        8: 2.309, 9: 2.414, 10: 2.458, 11: 2.581, 12: 2.732, 13: 2.828,
        14: 2.873, 15: 3.0, 16: 3.055, 17: 3.155, 18: 3.236, 19: 3.314,
        20: 3.314, 21: 3.464, 22: 3.5, 23: 3.581, 24: 3.618, 25: 3.732,
        26: 3.828, 27: 3.864, 28: 3.911, 29: 4.0, 30: 4.0, 31: 4.055,
        32: 4.155, 33: 4.215, 34: 4.236, 35: 4.309, 36: 4.328, 37: 4.414,
        38: 4.464, 39: 4.478, 40: 4.478
    }

    min_cyl_area = float('inf')
    best_cyl_config = None
    best_cyl_dims = None

    for n_layers in range(1, int(NUM_BALLS) + 1):
        n_per_layer = math.ceil(NUM_BALLS / n_layers)
        if n_per_layer > max(packing_ratios.keys()):
            continue

        ratio = packing_ratios[n_per_layer]
        
        # Calculate ideal radius and then round up to precision
        R_ideal = ratio * BALL_RADIUS
        R_cyl = math.ceil(R_ideal / PRECISION) * PRECISION
        
        H_cyl = n_layers * BALL_DIAMETER
        
        area = 2 * math.pi * R_cyl * (R_cyl + H_cyl)
        
        if area < min_cyl_area:
            min_cyl_area = area
            best_cyl_config = (n_layers, n_per_layer)
            best_cyl_dims = (R_cyl, H_cyl)

    cyl_material_cost = min_cyl_area * MATERIAL_COST_PER_CM2
    total_cyl_cost = cyl_material_cost + TOTAL_BALL_COST

    # --- Step 4: Compare and Output Final Answer ---
    print("--- Design Analysis ---")
    print(f"Minimum number of energy balls required: {int(NUM_BALLS)}")
    print(f"Fixed cost of energy balls: ${int(TOTAL_BALL_COST)}")
    print("-" * 25)
    
    if total_box_cost < total_cyl_cost:
        final_cost = total_box_cost
        material_cost = box_material_cost
        design_type = "Box"
        nx, ny, nz = best_box_config
        L, W, H = nx * BALL_DIAMETER, ny * BALL_DIAMETER, nz * BALL_DIAMETER
        print(f"Optimal Design: {design_type}")
        print(f"Configuration: {nx}x{ny}x{nz} balls")
        print(f"Dimensions (LxWxH): {L:.1f} cm x {W:.1f} cm x {H:.1f} cm")
        print(f"Surface Area: {min_box_area:.2f} cm^2")
    else:
        final_cost = total_cyl_cost
        material_cost = cyl_material_cost
        design_type = "Cylinder"
        n_layers, n_per_layer = best_cyl_config
        R_cyl, H_cyl = best_cyl_dims
        print(f"Optimal Design: {design_type}")
        print(f"Configuration: {n_layers} layers of up to {n_per_layer} balls")
        print(f"Dimensions (Radius x Height): {R_cyl:.1f} cm x {H_cyl:.1f} cm")
        print(f"Surface Area: {min_cyl_area:.2f} cm^2")

    print("-" * 25)
    print("--- Final Cost Calculation ---")
    print(f"Material Cost = Surface Area * Cost/cm^2")
    print(f"Material Cost = {min_cyl_area:.2f} cm^2 * ${MATERIAL_COST_PER_CM2}/cm^2 = ${material_cost:,.2f}")
    print(f"Total Cost = Material Cost + Ball Cost")
    print(f"Total Cost = ${material_cost:,.2f} + ${TOTAL_BALL_COST:,.2f}")
    print(f"Total Cost (C) = ${final_cost:,.2f}")
    
    # Return the final cost as an integer for the answer format
    return int(round(final_cost))

# --- Execute the solution ---
if __name__ == "__main__":
    final_answer = solve_container_problem()
    print(f"\n<<<C = {final_answer}>>>")
