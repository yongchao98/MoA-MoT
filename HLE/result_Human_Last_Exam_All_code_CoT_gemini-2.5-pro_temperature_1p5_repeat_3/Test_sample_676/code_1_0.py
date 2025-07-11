import math

def solve_container_problem():
    """
    Solves the container design optimization problem to find the minimum total cost.
    """

    # 1. Define Constants
    ENERGY_PER_BALL = 25  # MJ
    TARGET_ENERGY = 1000  # MJ
    BALL_COST = 1000      # USD
    MATERIAL_COST = 200   # USD/cm^2
    BALL_RADIUS = 2       # cm
    BALL_DIAMETER = 4     # cm
    PRECISION = 0.5       # cm

    # 2. Calculate Minimum Ball Count
    min_balls_needed = math.ceil(TARGET_ENERGY / ENERGY_PER_BALL)

    # --- Box Optimization ---
    min_cost_box = float('inf')
    best_box_config = {}

    # Iterate through possible numbers of balls to find a compact arrangement
    for n_actual in range(min_balls_needed, min_balls_needed * 2):
        # Find 3-integer factorizations of n_actual to get dimensions
        for nx in range(1, int(n_actual**(1/3.0)) + 2):
            if n_actual % nx == 0:
                for ny in range(nx, int(math.sqrt(n_actual / nx)) + 2):
                    if (n_actual / nx) % ny == 0:
                        nz = n_actual // (nx * ny)
                        
                        # Dimensions of the box
                        L, W, H = BALL_DIAMETER * nx, BALL_DIAMETER * ny, BALL_DIAMETER * nz
                        
                        # Calculate costs
                        sa = 2 * (L*W + L*H + W*H)
                        cost_balls = n_actual * BALL_COST
                        cost_material = sa * MATERIAL_COST
                        total_cost = cost_balls + cost_material
                        
                        if total_cost < min_cost_box:
                            min_cost_box = total_cost
                            best_box_config = {
                                "type": "Box",
                                "cost": total_cost,
                                "balls": n_actual,
                                "dims": (L, W, H),
                                "ball_dims": (nx,ny,nz),
                                "sa": sa,
                                "cost_balls": cost_balls,
                                "cost_material": cost_material
                            }

    # --- Cylinder Optimization ---
    min_cost_cylinder = float('inf')
    best_cylinder_config = {}

    # Using pre-computed optimal packing ratios for circles in a circle
    # C_k = Radius of large circle / radius of small circle
    PACKING_CONSTANTS = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.305, 9: 3.613, 10: 3.813, 11: 3.924, 12: 4.030, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 5.0,
        20: 5.122
    }

    for k, C_k in PACKING_CONSTANTS.items():
        # Number of layers needed
        m = math.ceil(min_balls_needed / k)
        n_actual = m * k
        
        # Calculate cylinder dimensions with precision
        H_cyl = m * BALL_DIAMETER
        R_inner = C_k * BALL_RADIUS
        R_cyl = math.ceil(R_inner / PRECISION) * PRECISION
        
        # Calculate costs
        sa = (2 * math.pi * R_cyl**2) + (2 * math.pi * R_cyl * H_cyl)
        cost_balls = n_actual * BALL_COST
        cost_material = sa * MATERIAL_COST
        total_cost = cost_balls + cost_material

        if total_cost < min_cost_cylinder:
            min_cost_cylinder = total_cost
            best_cylinder_config = {
                "type": "Cylinder",
                "cost": total_cost,
                "balls": n_actual,
                "dims": (R_cyl, H_cyl),
                "ball_dims": (k, m), # (balls per layer, layers)
                "sa": sa,
                "cost_balls": cost_balls,
                "cost_material": cost_material
            }
            
    # --- Compare and Conclude ---
    if min_cost_box < min_cost_cylinder:
        winner = best_box_config
    else:
        winner = best_cylinder_config

    print("Optimal Design Found:")
    print(f"Container Type: {winner['type']}")
    if winner['type'] == 'Box':
        L, W, H = winner['dims']
        nx, ny, nz = winner['ball_dims']
        print(f"Dimensions: Length={L} cm, Width={W} cm, Height={H} cm")
        print(f"Ball Arrangement: {nx} x {ny} x {nz}")
        print("\nCost Calculation:")
        print(f"Total Cost = (Balls * Cost/Ball) + (Surface Area * Cost/cm^2)")
        print(f"C = ({winner['balls']} * {BALL_COST}) + (2 * ({L}*{W} + {L}*{H} + {W}*{H})) * {MATERIAL_COST}")
    else:
        R, H = winner['dims']
        k, m = winner['ball_dims']
        print(f"Dimensions: Radius={R} cm, Height={H} cm")
        print(f"Ball Arrangement: {m} layers of {k} balls each")
        print("\nCost Calculation:")
        print(f"Total Cost = (Balls * Cost/Ball) + (Surface Area * Cost/cm^2)")
        print(f"C = ({winner['balls']} * {BALL_COST}) + (2 * pi * {R}^2 + 2 * pi * {R} * {H}) * {MATERIAL_COST}")

    final_cost = round(winner['cost'])
    print(f"\nLowest Total Cost (C): ${winner['cost']:.2f}")
    print(f"Final Answer (rounded to nearest dollar): C = {final_cost}")
    print(f"<<<{final_cost}>>>")


solve_container_problem()
