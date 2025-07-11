import math

def solve_container_problem():
    """
    Solves the container design optimization problem to find the minimum total cost.
    """
    # --- Constants ---
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = 4.0 # cm
    BALL_ENERGY = 25  # MJ
    BALL_COST = 1000  # usd
    MATERIAL_COST = 200  # usd per cm^2
    MIN_TOTAL_ENERGY = 1000  # MJ
    PRECISION = 0.5 # cm

    # --- Step 1: Calculate minimum number of balls ---
    min_balls = math.ceil(MIN_TOTAL_ENERGY / BALL_ENERGY)

    # --- Initialize variables to track the best design ---
    min_total_cost = float('inf')
    best_design_details = ""

    def get_factor_triplets(n):
        """Helper function to find unique, sorted integer factor triplets of n."""
        triplets = set()
        for i in range(1, int(n**0.5) + 2):
            if n % i == 0:
                n_div_i = n // i
                for j in range(i, int(n_div_i**0.5) + 2):
                    if n_div_i % j == 0:
                        k = n_div_i // j
                        if k >= j:
                            triplets.add(tuple(sorted((i, j, k))))
        if not triplets:
             # Handle prime numbers
             triplets.add((1,1,n))
        return list(triplets)

    # --- Search through a range of ball numbers ---
    # The cost of balls increases linearly, so the optimum is likely near the minimum.
    # We search a bit beyond the minimum just in case a better packing reduces surface area significantly.
    for n_balls in range(min_balls, min_balls + 40):
        ball_cost = n_balls * BALL_COST
        factor_triplets = get_factor_triplets(n_balls)

        # --- Box Design Optimization ---
        for nx, ny, nz in factor_triplets:
            l, w, h = nx * BALL_DIAMETER, ny * BALL_DIAMETER, nz * BALL_DIAMETER
            surface_area = 2 * (l*w + l*h + w*h)
            container_cost = surface_area * MATERIAL_COST
            total_cost = ball_cost + container_cost

            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design_details = (
                    f"The optimal design is a Box container.\n\n"
                    f"Equation for the total cost C:\n"
                    f"C = (Number of balls * Cost per ball) + (Container surface area * Material cost)\n\n"
                    f"Calculation:\n"
                    f"Number of balls (N) = {n_balls}\n"
                    f"Arrangement of balls = {nx} x {ny} x {nz}\n"
                    f"Cost of balls = {n_balls} * {BALL_COST} = {ball_cost}\n\n"
                    f"Container length (L) = {nx} * {BALL_DIAMETER} = {l:.1f} cm\n"
                    f"Container width (W) = {ny} * {BALL_DIAMETER} = {w:.1f} cm\n"
                    f"Container height (H) = {nz} * {BALL_DIAMETER} = {h:.1f} cm\n"
                    f"Container surface area = 2 * (L*W + L*H + W*H)\n"
                    f"                     = 2 * ({l:.1f}*{w:.1f} + {l:.1f}*{h:.1f} + {w:.1f}*{h:.1f}) = {surface_area:.2f} cm^2\n"
                    f"Cost of container = {surface_area:.2f} * {MATERIAL_COST} = {container_cost:.2f}\n\n"
                    f"Total Cost (C) = {ball_cost} + {container_cost:.2f} = {total_cost:.2f}"
                )

        # --- Cylinder Design Optimization ---
        for factors in factor_triplets:
            # Test each factor as the number of layers 'm'
            for i in range(3):
                m = factors[i]
                nx, ny = sorted((factors[(i + 1) % 3], factors[(i + 2) % 3]))

                h_cyl = m * BALL_DIAMETER
                
                # Required radius for a grid of spheres
                if nx == 1 and ny == 1:
                    r_req = BALL_RADIUS
                else:
                    dist_to_center = 0.5 * BALL_DIAMETER * math.sqrt((nx - 1)**2 + (ny - 1)**2)
                    r_req = dist_to_center + BALL_RADIUS

                # Apply precision constraint, rounding up to the nearest 0.5 cm
                r_cyl = math.ceil(r_req / PRECISION) * PRECISION

                surface_area_cyl = 2 * math.pi * r_cyl * (r_cyl + h_cyl)
                container_cost_cyl = surface_area_cyl * MATERIAL_COST
                total_cost_cyl = ball_cost + container_cost_cyl

                if total_cost_cyl < min_total_cost:
                    min_total_cost = total_cost_cyl
                    best_design_details = (
                        f"The optimal design is a Cylinder container.\n\n"
                        f"Equation for the total cost C:\n"
                        f"C = (Number of balls * Cost per ball) + (Container surface area * Material cost)\n\n"
                        f"Calculation:\n"
                        f"Number of balls (N) = {n_balls}\n"
                        f"Arrangement of balls = {m} layers of {nx} x {ny} grids\n"
                        f"Cost of balls = {n_balls} * {BALL_COST} = {ball_cost}\n\n"
                        f"Container height (H) = {m} * {BALL_DIAMETER} = {h_cyl:.1f} cm\n"
                        f"Container radius (R) = {r_cyl:.1f} cm (rounded up for precision)\n"
                        f"Container surface area = 2 * pi * R * (R + H)\n"
                        f"                     = 2 * pi * {r_cyl:.1f} * ({r_cyl:.1f} + {h_cyl:.1f}) = {surface_area_cyl:.2f} cm^2\n"
                        f"Cost of container = {surface_area_cyl:.2f} * {MATERIAL_COST} = {container_cost_cyl:.2f}\n\n"
                        f"Total Cost (C) = {ball_cost} + {container_cost_cyl:.2f} = {total_cost_cyl:.2f}"
                    )

    if min_total_cost == float('inf'):
        print("Could not find a valid design.")
        print("<<<0>>>")
    else:
        print(best_design_details)
        print(f"\n<<<{min_total_cost:.0f}>>>")

solve_container_problem()