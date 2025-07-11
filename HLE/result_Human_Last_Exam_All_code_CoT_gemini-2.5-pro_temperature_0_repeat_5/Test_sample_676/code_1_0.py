import math

def solve_pioneer_probe_design():
    """
    Calculates the minimum cost to design a container for energy balls,
    considering both box and cylinder shapes.
    """

    # --- Problem Constants ---
    ENERGY_REQ_MJ = 1000
    ENERGY_PER_BALL_MJ = 25
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = BALL_RADIUS_CM * 2
    BALL_COST_USD = 1000
    MATERIAL_COST_PER_CM2_USD = 200
    PRECISION_CM = 0.5

    # --- Step 1: Calculate minimum number of balls ---
    num_balls_req = math.ceil(ENERGY_REQ_MJ / ENERGY_PER_BALL_MJ)

    # --- Step 2: Analyze Box Container ---
    min_box_cost = float('inf')
    best_box_config = {}

    # To minimize surface area for a given volume, dimensions should be as close as possible.
    # We test integer factorizations of 40.
    factorizations = [
        (40, 1, 1), (20, 2, 1), (10, 4, 1), (10, 2, 2), (8, 5, 1), (5, 4, 2)
    ]

    for nx, ny, nz in factorizations:
        num_balls = nx * ny * nz
        if num_balls >= num_balls_req:
            L = nx * BALL_DIAMETER_CM
            W = ny * BALL_DIAMETER_CM
            H = nz * BALL_DIAMETER_CM
            
            area = 2 * (L*W + W*H + H*L)
            material_cost = area * MATERIAL_COST_PER_CM2_USD
            ball_cost = num_balls * BALL_COST_USD
            total_cost = material_cost + ball_cost

            if total_cost < min_box_cost:
                min_box_cost = total_cost
                best_box_config = {
                    "type": "Box",
                    "dims": (L, W, H),
                    "balls": num_balls,
                    "area": area,
                    "cost": total_cost
                }

    # --- Step 3: Analyze Cylinder Container ---
    min_cyl_cost = float('inf')
    best_cyl_config = {}

    # Packing ratios (R_container / r_ball) for packing n circles in a circle.
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cinc/cinc.html
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0, 8: 3.304, 9: 3.613,
        10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236, 14: 4.328, 15: 4.521, 16: 4.615,
        17: 4.792, 18: 4.864, 19: 4.864, 20: 5.122, 21: 5.248, 22: 5.323, 23: 5.532,
        24: 5.607, 25: 5.685, 26: 5.864, 27: 5.889, 28: 6.092, 29: 6.122, 30: 6.246,
        31: 6.398, 32: 6.468, 33: 6.519, 34: 6.691, 35: 6.747, 36: 6.747, 37: 6.747,
        38: 6.999, 39: 7.143, 40: 7.185
    }

    for n_cols in range(1, num_balls_req + 1):
        n_rows = math.ceil(num_balls_req / n_cols)
        actual_balls = n_cols * n_rows
        
        H = n_rows * BALL_DIAMETER_CM
        
        # Calculate required radius for the base
        ratio = packing_ratios[n_cols]
        R_pack = BALL_RADIUS_CM * ratio
        
        # Apply manufacturing precision
        R = math.ceil(R_pack / PRECISION_CM) * PRECISION_CM
        
        area = 2 * math.pi * R * H + 2 * math.pi * R**2
        material_cost = area * MATERIAL_COST_PER_CM2_USD
        ball_cost = actual_balls * BALL_COST_USD
        total_cost = material_cost + ball_cost
        
        if total_cost < min_cyl_cost:
            min_cyl_cost = total_cost
            best_cyl_config = {
                "type": "Cylinder",
                "dims": (R, H),
                "balls": actual_balls,
                "area": area,
                "cost": total_cost,
                "cols": n_cols
            }

    # --- Step 4: Compare and Conclude ---
    if min_box_cost < min_cyl_cost:
        best_design = best_box_config
        L, W, H = best_design['dims']
        print(f"The optimal container is a Box with dimensions L={L:.1f} cm, W={W:.1f} cm, H={H:.1f} cm.")
        print(f"This container holds {best_design['balls']} balls.")
        print(f"The surface area is 2 * ({L:.1f}*{W:.1f} + {W:.1f}*{H:.1f} + {H:.1f}*{L:.1f}) = {best_design['area']:.2f} cm^2.")
        material_cost = best_design['area'] * MATERIAL_COST_PER_CM2_USD
        ball_cost = best_design['balls'] * BALL_COST_USD
        print(f"The total cost is ({best_design['area']:.2f} cm^2 * ${MATERIAL_COST_PER_CM2_USD}) + ({best_design['balls']} balls * ${BALL_COST_USD}) = ${material_cost:.0f} + ${ball_cost:.0f} = ${best_design['cost']:.0f}.")
        final_cost = round(best_design['cost'])
    else:
        best_design = best_cyl_config
        R, H = best_design['dims']
        print(f"The optimal container is a Cylinder with radius R={R:.1f} cm and height H={H:.1f} cm.")
        print(f"This container holds {best_design['balls']} balls (arranged in {best_design['cols']} columns).")
        print(f"The surface area is 2*pi*{R:.1f}^2 + 2*pi*{R:.1f}*{H:.1f} = {best_design['area']:.2f} cm^2.")
        material_cost = best_design['area'] * MATERIAL_COST_PER_CM2_USD
        ball_cost = best_design['balls'] * BALL_COST_USD
        print(f"The total cost is ({best_design['area']:.2f} cm^2 * ${MATERIAL_COST_PER_CM2_USD}) + ({best_design['balls']} balls * ${BALL_COST_USD}) = ${material_cost:.0f} + ${ball_cost:.0f} = ${round(best_design['cost']):.0f}.")
        final_cost = round(best_design['cost'])

    print(f"\nFinal Answer: The lowest total cost is ${final_cost}.")
    print(f"<<<{final_cost}>>>")

solve_pioneer_probe_design()