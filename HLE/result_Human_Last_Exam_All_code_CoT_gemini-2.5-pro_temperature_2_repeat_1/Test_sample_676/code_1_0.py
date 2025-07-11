import math

def solve_pioneer_probe_design():
    """
    Calculates the minimum cost to design a container for energy balls,
    considering both box and cylinder shapes.
    """
    # 1. Define constants based on the problem description
    MIN_ENERGY_MJ = 1000
    ENERGY_PER_BALL_MJ = 25
    BALL_COST_USD = 1000
    BALL_DIAMETER_CM = 4
    MATERIAL_COST_PER_CM2_USD = 200
    PRECISION_CM = 0.5

    # 2. Determine the minimum number of energy balls
    # The total cost function increases with the number of balls, so the optimum
    # will be at the minimum required number.
    n_balls = math.ceil(MIN_ENERGY_MJ / ENERGY_PER_BALL_MJ)

    min_total_cost = float('inf')
    best_config = {}

    # 3. Find all unique integer factor triplets for the number of balls
    def get_factor_triplets(n):
        triplets = set()
        # To find factors a, b, c of n: a <= b <= c
        # We must have a <= n^(1/3)
        for i in range(1, int(n**(1/3.0)) + 2):
            if n % i == 0:
                n_div_i = n // i
                # And b <= sqrt(n/a)
                for j in range(i, int(n_div_i**0.5) + 2):
                    if n_div_i % j == 0:
                        k = n_div_i // j
                        # Sort to handle duplicates (e.g., (2,4,5) is same as (5,2,4))
                        triplets.add(tuple(sorted((i, j, k))))
        return [t for t in triplets if t[0] * t[1] * t[2] == n]

    factor_triplets = get_factor_triplets(n_balls)

    # 4. Iterate through each arrangement to find the optimal container
    for triplet in factor_triplets:
        nx, ny, nz = triplet

        # --- Box Container Calculation ---
        L = nx * BALL_DIAMETER_CM
        W = ny * BALL_DIAMETER_CM
        H = nz * BALL_DIAMETER_CM
        
        sa_box = 2 * (L * W + L * H + W * H)
        cost_balls = n_balls * BALL_COST_USD
        cost_material_box = sa_box * MATERIAL_COST_PER_CM2_USD
        total_cost_box = cost_balls + cost_material_box
        
        if total_cost_box < min_total_cost:
            min_total_cost = total_cost_box
            best_config = {
                "type": "Box", "n_balls": n_balls, "dims_balls": (nx, ny, nz),
                "dims_cm": (L, W, H), "surface_area": sa_box,
                "cost_balls": cost_balls, "cost_material": cost_material_box,
                "total_cost": total_cost_box
            }

        # --- Cylinder Container Calculation (check all 3 orientations) ---
        import itertools
        # Use set to avoid re-calculating for symmetric triplets like (2,2,10)
        for p in set(itertools.permutations(triplet)):
            base_dim1, base_dim2, h_dim = p[0], p[1], p[2]
            
            H_cyl = h_dim * BALL_DIAMETER_CM
            
            # The rectangular base of balls (base_dim1 x base_dim2) must fit in a circle
            # The diameter of that circle is the diagonal of the rectangle
            required_radius = 0.5 * BALL_DIAMETER_CM * math.sqrt(base_dim1**2 + base_dim2**2)
            
            # The container's radius must adhere to the manufacturing precision
            container_radius = math.ceil(required_radius / PRECISION_CM) * PRECISION_CM
            
            sa_cyl = 2 * math.pi * container_radius**2 + 2 * math.pi * container_radius * H_cyl
            cost_material_cyl = sa_cyl * MATERIAL_COST_PER_CM2_USD
            total_cost_cyl = cost_balls + cost_material_cyl

            if total_cost_cyl < min_total_cost:
                min_total_cost = total_cost_cyl
                best_config = {
                    "type": "Cylinder", "n_balls": n_balls, "dims_balls": (base_dim1, base_dim2, h_dim),
                    "dims_cm": (container_radius, H_cyl), "surface_area": sa_cyl,
                    "cost_balls": cost_balls, "cost_material": cost_material_cyl,
                    "total_cost": total_cost_cyl
                }
    
    # 5. Output the final result and the equation
    if not best_config:
        print("A solution could not be found.")
        print("<<<0>>>")
        return

    n_final = best_config["n_balls"]
    cost_balls_final = best_config["cost_balls"]
    sa_final = best_config["surface_area"]
    cost_material_final = best_config["cost_material"]
    cost_total_final = best_config["total_cost"]

    print("--- Optimal Design ---")
    print(f"The best design is a {best_config['type']} container.")
    if best_config['type'] == 'Box':
        l, w, h = best_config['dims_cm']
        print(f"Dimensions: {l:.1f} cm x {w:.1f} cm x {h:.1f} cm")
    else:
        r, h = best_config['dims_cm']
        print(f"Dimensions: Radius {r:.1f} cm, Height {h:.1f} cm")
    
    print("\n--- Final Cost Calculation ---")
    print("Total Cost = (Number of Balls * Cost per Ball) + (Surface Area * Cost per cm^2)")
    print(f"Total Cost = ({n_final} * {BALL_COST_USD}) + ({sa_final:.2f} * {MATERIAL_COST_PER_CM2_USD})")
    print(f"Total Cost = {int(cost_balls_final)} + {int(round(cost_material_final))} = {int(round(cost_total_final))}")

    print(f"\n<<<{int(round(cost_total_final))}>>>")

solve_pioneer_probe_design()