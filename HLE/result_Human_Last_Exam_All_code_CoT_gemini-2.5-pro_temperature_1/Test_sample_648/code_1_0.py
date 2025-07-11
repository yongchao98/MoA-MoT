import math

def solve_pioneer_probe_design():
    """
    Solves the container design optimization problem for the Pioneer probe.
    It searches for the minimum cost box or cylinder container that can hold
    at least 1000 MJ of energy, subject to material and manufacturing constraints.
    """
    # --- Constants ---
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = 4.0  # cm
    BALL_ENERGY = 30  # MJ
    BALL_COST = 1000  # USD
    MATERIAL_COST_PER_CM2 = 200  # USD
    MAX_SURFACE_AREA = 1000.0  # cm^2
    ENERGY_REQUIREMENT = 1000.0  # MJ
    PRECISION = 0.5  # cm

    # Step 1: Calculate minimum number of balls
    min_balls_needed = math.ceil(ENERGY_REQUIREMENT / BALL_ENERGY)

    # --- Best solution tracking ---
    best_cost = float('inf')
    best_design = None

    # --- Part 1: Box container search ---
    # l, w, h are in units of PRECISION (0.5 cm)
    # To fit at least one ball, L,W,H must be >= 4cm, so l,w,h >= 4/0.5 = 8
    min_dim_units = int(BALL_DIAMETER / PRECISION)
    # Set a reasonable search limit for dimensions (e.g., up to 25 cm)
    search_limit_units = int(25 / PRECISION)

    for l in range(min_dim_units, search_limit_units + 1):
        for w in range(min_dim_units, l + 1):
            for h in range(min_dim_units, w + 1):
                L, W, H = l * PRECISION, w * PRECISION, h * PRECISION
                surface_area = 2 * (L * W + W * H + H * L)

                if surface_area > MAX_SURFACE_AREA:
                    if h == min_dim_units: # Optimization: if the smallest h fails, larger w will also fail
                        if w == min_dim_units: l_will_fail = True
                        else: w_will_fail = True
                    break # Since h is increasing, any larger h will also exceed SA
                
                num_balls = (math.floor(L / BALL_DIAMETER) *
                             math.floor(W / BALL_DIAMETER) *
                             math.floor(H / BALL_DIAMETER))

                if num_balls >= min_balls_needed:
                    total_cost = (surface_area * MATERIAL_COST_PER_CM2) + (num_balls * BALL_COST)
                    if total_cost < best_cost:
                        best_cost = total_cost
                        best_design = {
                            "type": "Box", "L": L, "W": W, "H": H,
                            "SA": surface_area, "N_balls": num_balls, "cost": total_cost
                        }
            if 'w_will_fail' in locals() and w_will_fail: break
        if 'l_will_fail' in locals() and l_will_fail: break


    # --- Part 2: Cylinder container search ---
    # Data for optimal packing of k circles in a larger circle: k -> R_container/r_ball
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.502, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 5.0, 20: 5.122
    }

    for k in sorted(packing_ratios.keys()):
        num_layers = math.ceil(min_balls_needed / k)
        H = num_layers * BALL_DIAMETER

        min_R = BALL_RADIUS * packing_ratios[k]
        R = math.ceil(min_R / PRECISION) * PRECISION

        surface_area = 2 * math.pi * R**2 + 2 * math.pi * R * H

        if surface_area <= MAX_SURFACE_AREA:
            num_balls = num_layers * k
            total_cost = (surface_area * MATERIAL_COST_PER_CM2) + (num_balls * BALL_COST)
            if total_cost < best_cost:
                best_cost = total_cost
                best_design = {
                    "type": "Cylinder", "R": R, "H": H,
                    "SA": surface_area, "N_balls": num_balls, "cost": total_cost
                }

    # --- Part 3: Output result ---
    if best_design is None:
        final_cost = 0
        print("No valid container design found that meets all constraints.")
    else:
        final_cost = best_design['cost']
        print(f"Optimal Design Found: {best_design['type']} Container")
        print("-" * 30)
        if best_design['type'] == 'Box':
            print(f"Dimensions: Length = {best_design['L']} cm, Width = {best_design['W']} cm, Height = {best_design['H']} cm")
        else:  # Cylinder
            print(f"Dimensions: Radius = {best_design['R']} cm, Height = {best_design['H']} cm")
        
        print(f"Container Surface Area: {best_design['SA']:.2f} cm^2")
        print(f"Energy Balls Packed: {best_design['N_balls']}")
        print(f"Total Energy Provided: {best_design['N_balls'] * BALL_ENERGY} MJ (Requirement: >= {ENERGY_REQUIREMENT} MJ)")
        print("-" * 30)
        print("Final Cost Calculation:")
        
        cost_material = best_design['SA'] * MATERIAL_COST_PER_CM2
        cost_balls = best_design['N_balls'] * BALL_COST
        
        print(f"C = (Surface Area * Material Cost) + (Number of Balls * Ball Cost)")
        print(f"C = ({best_design['SA']:.2f} cm^2 * {MATERIAL_COST_PER_CM2} usd/cm^2) + ({best_design['N_balls']} balls * {BALL_COST} usd/ball)")
        print(f"C = {cost_material:.2f} + {cost_balls:.2f}")
        print(f"C = {best_design['cost']:.2f} usd")

    print(f"\n<<<{final_cost:.2f}>>>")

# Execute the solution
solve_pioneer_probe_design()