import math

def solve_pioneer_probe_design():
    """
    Solves the container design optimization problem for the Pioneer probe.
    It searches for the lowest cost configuration (box or cylinder) that can hold
    at least 1000 MJ of energy, with a surface area of at most 1000 cm^2.
    """
    # --- Problem Parameters ---
    MIN_ENERGY = 1000          # MJ
    ENERGY_PER_BALL = 30       # MJ
    BALL_COST = 1000           # USD per ball
    MAX_SA = 1000              # cm^2
    MATERIAL_COST = 200        # USD per cm^2
    PRECISION = 0.5            # cm
    BALL_DIAMETER = 4          # cm

    min_balls_needed = math.ceil(MIN_ENERGY / ENERGY_PER_BALL)

    min_total_cost = float('inf')
    best_design = None

    # --- 1. Box Container Search ---
    # Iterate through possible dimensions (L, W, H) as multiples of PRECISION.
    # We use loops with optimized ranges based on geometric constraints.
    # L_max is ~sqrt(1000/2) = 22.3 cm. integer_l max is ~45.
    l_min_int = math.ceil(BALL_DIAMETER / PRECISION)
    l_max_int = math.ceil(math.sqrt(MAX_SA / 2) / PRECISION)

    for l_int in range(l_min_int, l_max_int + 2):
        L = l_int * PRECISION
        for w_int in range(l_int, math.ceil((MAX_SA / (2 * L)) / PRECISION) + 1):
            W = w_int * PRECISION
            if 2 * L * W > MAX_SA:
                continue
            
            h_max_val = (MAX_SA / 2 - L * W) / (L + W) if (L + W) > 0 else 0
            if h_max_val < W:
                continue
            h_max_int = math.floor(h_max_val / PRECISION)
            
            for h_int in range(w_int, h_max_int + 1):
                H = h_int * PRECISION
                surface_area = 2 * (L*W + W*H + H*L)
                if surface_area > MAX_SA:
                    continue

                total_balls = math.floor(L / BALL_DIAMETER) * math.floor(W / BALL_DIAMETER) * math.floor(H / BALL_DIAMETER)

                if total_balls >= min_balls_needed:
                    current_cost = (surface_area * MATERIAL_COST) + (total_balls * BALL_COST)
                    if current_cost < min_total_cost:
                        min_total_cost = current_cost
                        best_design = {"type": "Box", "L": L, "W": W, "H": H,
                                       "surface_area": surface_area, "num_balls": total_balls, "cost": current_cost}

    # --- 2. Cylinder Container Search ---
    # Iterate through possible dimensions (R, H) as multiples of PRECISION.
    # r_max is ~sqrt(1000 / (2*pi)) = 12.6 cm. integer_r max is ~25.
    r_max_int = math.ceil(math.sqrt(MAX_SA / (2 * math.pi)) / PRECISION)
    r_min_int = math.ceil((BALL_DIAMETER / math.sqrt(2)) / PRECISION)

    for r_int in range(r_min_int, r_max_int + 1):
        R = r_int * PRECISION
        h_max_val = (MAX_SA / (2 * math.pi * R)) - R if R > 0 else 0
        if h_max_val < BALL_DIAMETER:
            continue
        h_max_int = math.floor(h_max_val / PRECISION)
        
        for h_int in range(l_min_int, h_max_int + 1):
            H = h_int * PRECISION
            surface_area = 2 * math.pi * R * (R + H)
            if surface_area > MAX_SA:
                continue
            
            num_balls_per_layer = math.floor(math.sqrt(2) * R / BALL_DIAMETER)**2
            num_layers = math.floor(H / BALL_DIAMETER)
            total_balls = num_balls_per_layer * num_layers

            if total_balls >= min_balls_needed:
                current_cost = (surface_area * MATERIAL_COST) + (total_balls * BALL_COST)
                if current_cost < min_total_cost:
                    min_total_cost = current_cost
                    best_design = {"type": "Cylinder", "R": R, "H": H,
                                   "surface_area": surface_area, "num_balls": total_balls, "cost": current_cost}

    # --- 3. Final Output ---
    if best_design is None:
        final_cost = 0
        print("No feasible design found.")
    else:
        final_cost = round(best_design['cost'])
        container_cost = best_design['surface_area'] * MATERIAL_COST
        balls_cost = best_design['num_balls'] * BALL_COST
        # Output the numbers in the final equation as requested
        print(f"Total Cost = {container_cost:.0f} + {balls_cost:.0f} = {final_cost}")
        
    print(f"<<<{final_cost}>>>")

solve_pioneer_probe_design()