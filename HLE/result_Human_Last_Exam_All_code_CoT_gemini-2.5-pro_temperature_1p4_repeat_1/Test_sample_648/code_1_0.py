import math

def solve_container_design():
    # --- Problem Parameters ---
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = BALL_RADIUS * 2
    BALL_ENERGY = 30  # MJ
    BALL_COST = 1000  # USD
    MAX_SA = 1000.0  # cm^2
    MATERIAL_COST_PER_CM2 = 200  # USD
    MIN_TOTAL_ENERGY = 1000  # MJ
    PRECISION = 0.5  # cm

    # --- Derived Parameters ---
    MIN_REQUIRED_BALLS = math.ceil(MIN_TOTAL_ENERGY / BALL_ENERGY)

    min_total_cost = float('inf')
    best_design = None

    # --- Search for the best Box design ---
    # The search range is determined by the max SA. If L=W=H, 6L^2<=1000 -> L<=12.9.
    # So dimension units (in 0.5cm) go up to ~26. We'll search a bit more.
    max_dim_units = 35 
    for h_units in range(1, max_dim_units + 1):
        H = h_units * PRECISION
        for w_units in range(1, h_units + 1):
            W = w_units * PRECISION
            for l_units in range(1, w_units + 1):
                L = l_units * PRECISION

                sa = 2 * (L * W + L * H + W * H)
                if sa > MAX_SA:
                    continue

                n_balls = (math.floor(L / BALL_DIAMETER) *
                           math.floor(W / BALL_DIAMETER) *
                           math.floor(H / BALL_DIAMETER))

                if n_balls < MIN_REQUIRED_BALLS:
                    continue

                material_cost = sa * MATERIAL_COST_PER_CM2
                ball_cost = n_balls * BALL_COST
                total_cost = material_cost + ball_cost

                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        "type": "Box",
                        "dims": {"L": L, "W": W, "H": H},
                        "n_balls": n_balls,
                        "sa": sa,
                        "cost": total_cost,
                        "ball_cost": ball_cost,
                        "material_cost": material_cost
                    }

    # --- Search for the best Cylinder design ---
    # Search range for radius (R): 2*pi*R^2 <= 1000 -> R <= 12.6.
    # Radius units up to 12.5/0.5 = 25. Let's search up to 30.
    max_r_units = 30
    for r_units in range(1, max_r_units + 1):
        R = r_units * PRECISION
        
        # Max height H for a given R to keep SA <= 1000
        # 2*pi*R^2 + 2*pi*R*H <= 1000 => H <= (1000 - 2*pi*R^2)/(2*pi*R)
        if 2 * math.pi * R**2 > MAX_SA:
            continue
        max_h = (MAX_SA - 2 * math.pi * R**2) / (2 * math.pi * R)
        max_h_units = math.floor(max_h / PRECISION)

        for h_units in range(1, max_h_units + 1):
            H = h_units * PRECISION
            
            sa = 2 * math.pi * R**2 + 2 * math.pi * R * H
            
            # Using area-ratio packing model
            n_layers = math.floor(H / BALL_DIAMETER)
            if n_layers == 0:
                continue
            
            n_per_layer = math.floor((math.pi * R**2) / (BALL_DIAMETER**2))
            if n_per_layer == 0:
                continue

            n_balls = n_layers * n_per_layer

            if n_balls < MIN_REQUIRED_BALLS:
                continue

            material_cost = sa * MATERIAL_COST_PER_CM2
            ball_cost = n_balls * BALL_COST
            total_cost = material_cost + ball_cost
            
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design = {
                    "type": "Cylinder",
                    "dims": {"R": R, "H": H},
                    "n_balls": n_balls,
                    "sa": sa,
                    "cost": total_cost,
                    "ball_cost": ball_cost,
                    "material_cost": material_cost
                }

    # --- Output the result ---
    if best_design:
        print("Optimal Design Found:")
        print(f"  Container Type: {best_design['type']}")
        if best_design['type'] == 'Box':
            dims = best_design['dims']
            print(f"  Dimensions: Length={dims['L']:.1f} cm, Width={dims['W']:.1f} cm, Height={dims['H']:.1f} cm")
        else:
            dims = best_design['dims']
            print(f"  Dimensions: Radius={dims['R']:.1f} cm, Height={dims['H']:.1f} cm")
        
        print(f"  Number of Energy Balls: {best_design['n_balls']}")
        print(f"  Total Energy: {best_design['n_balls'] * BALL_ENERGY} MJ")
        print(f"  Container Surface Area: {best_design['sa']:.2f} cm^2")
        print("\nCost Breakdown:")
        print(f"  Cost of Energy Balls: {best_design['n_balls']} balls * ${BALL_COST}/ball = ${best_design['ball_cost']:.2f}")
        print(f"  Cost of Material: {best_design['sa']:.2f} cm^2 * ${MATERIAL_COST_PER_CM2}/cm^2 = ${best_design['material_cost']:.2f}")
        print("\nFinal Equation for Total Cost (C):")
        print(f"C = ({best_design['n_balls']} * {BALL_COST}) + ({best_design['sa']:.2f} * {MATERIAL_COST_PER_CM2})")
        print(f"C = ${best_design['cost']:.2f}")
        
        # Storing final answer for the required format
        global final_answer_c
        final_answer_c = best_design['cost']
    else:
        print("No solution found that meets all criteria.")
        final_answer_c = 0

# Run the solver and print the final answer in the required format
final_answer_c = 0
solve_container_design()
print(f"\n<<<C = {final_answer_c:.2f}>>>")