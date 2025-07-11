import math

def solve_box_optimization():
    """
    This function seeks to find a more material-efficient container for energy balls.

    It works by first establishing the performance of the current container (a 12x12x12 cm cube)
    in terms of the number of balls it can hold and its surface area.

    Then, it formulates the problem as finding a new set of dimensions that can hold at
    least the same number of balls but with a smaller surface area. The problem is simplified
    by assuming an efficient simple cubic lattice packing for the spherical balls, which is
    well-suited for the grid-based constraints (all measurements are multiples of 0.5 cm).

    The core of the program is a search for integer ball counts along each axis (nx, ny, nz)
    that satisfy two key inequalities derived from the problem's constraints:
    1. The total number of balls (nx * ny * nz) must be at least that of the original box.
    2. The new surface area, which is a function of nx, ny, and nz, must be less than
       the original box's surface area.

    The program iterates through possible values for (nx, ny, nz) and checks if any
    combination meets these criteria. If a better box is found, its dimensions and surface
    area are printed. Otherwise, it concludes that no such box exists and prints '0'.
    """

    # --- Step 1: Analyze the initial container ---
    initial_dims = (12, 12, 12)
    ball_diameter = 4.0

    # Calculate the number of balls that fit in the initial container
    # using simple lattice packing.
    nx_initial = math.floor(initial_dims[0] / ball_diameter)
    ny_initial = math.floor(initial_dims[1] / ball_diameter)
    nz_initial = math.floor(initial_dims[2] / ball_diameter)
    min_balls_to_pack = nx_initial * ny_initial * nz_initial

    # Calculate the surface area of the initial container.
    l, w, h = initial_dims
    initial_surface_area = 2 * (l * w + l * h + w * h)

    # --- Step 2: Search for a better container ---
    # The problem reduces to finding integers nx, ny, nz such that:
    # 1. nx * ny * nz >= min_balls_to_pack (which is 27)
    # 2. 32 * (nx*ny + nx*nz + ny*nz) < initial_surface_area (which is 864)
    #    This simplifies to: nx*ny + nx*nz + ny*nz < 27

    best_solution = None

    # We search for triplets (nx, ny, nz). The search range can be small because
    # the second condition (nx*ny + nx*nz + ny*nz < 27) limits how large
    # nx, ny, and nz can be. A limit of 30 is very safe.
    search_limit = 30
    for nx in range(1, search_limit):
        # To avoid permutations and redundant checks, we enforce nx <= ny <= nz
        for ny in range(nx, search_limit):
            # Optimization: If the smallest possible sum already violates the condition,
            # we can break the inner loops early.
            if (nx * ny + nx * ny + ny * ny) >= 27:
                break
            for nz in range(ny, search_limit):
                # Condition 1: Check if it holds enough balls
                if (nx * ny * nz) < min_balls_to_pack:
                    continue

                # Condition 2: Check if the surface area is smaller
                pairwise_product_sum = nx * ny + nx * nz + ny * nz
                if pairwise_product_sum >= 27:
                    # Since nz is increasing, all further values in this loop will also fail.
                    break

                # If we reach here, a better solution is theoretically possible.
                # In this specific problem, this block will not be reached.
                l_new = nx * ball_diameter
                w_new = ny * ball_diameter
                h_new = nz * ball_diameter
                current_surface_area = 32 * pairwise_product_sum

                if best_solution is None or current_surface_area < best_solution[3]:
                    best_solution = (int(l_new), int(w_new), int(h_new), int(current_surface_area))


    # --- Step 3: Output the result ---
    if best_solution:
        l_final, w_final, h_final, sa_final = best_solution
        print(f"{l_final}:{w_final}:{h_final}:{sa_final}")
    else:
        print("0")

solve_box_optimization()