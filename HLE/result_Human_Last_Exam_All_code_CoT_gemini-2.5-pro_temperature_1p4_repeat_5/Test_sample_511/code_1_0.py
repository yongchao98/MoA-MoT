import math

def solve_packing_problem():
    """
    Finds an optimal container for energy balls, minimizing surface area
    while holding at least as many balls as the original container.
    """

    # --- Problem Configuration ---
    ORIGINAL_DIM = 12.0
    BALL_RADIUS = 2.0
    PRECISION = 0.5

    # --- Derived Constants ---
    BALL_DIAMETER_SQ = (2 * BALL_RADIUS)**2
    ORIGINAL_SURFACE_AREA = 6 * ORIGINAL_DIM**2

    # Caching to speed up repeated calls with the same dimensions
    packer_cache = {}

    def greedy_packer(l, w, h):
        """
        Calculates the number of non-overlapping spheres that can be packed
        into a box of size l x w x h using a greedy algorithm.
        """
        # Sort dimensions to ensure cache hits for permutations like 10x12x14 and 14x10x12
        dims = tuple(sorted((l, w, h)))
        if dims in packer_cache:
            return packer_cache[dims]

        min_coord = BALL_RADIUS
        max_x, max_y, max_z = dims[0] - BALL_RADIUS, dims[1] - BALL_RADIUS, dims[2] - BALL_RADIUS

        if dims[0] < 2 * BALL_RADIUS or dims[1] < 2 * BALL_RADIUS or dims[2] < 2 * BALL_RADIUS:
            return 0

        # Generate potential center points based on the 0.5cm precision requirement.
        # The list is implicitly sorted by z, then y, then x, which defines the greedy order.
        potential_centers = []
        z = min_coord
        while z <= max_z:
            y = min_coord
            while y <= max_y:
                x = min_coord
                while x <= max_x:
                    potential_centers.append((x, y, z))
                    x += PRECISION
                y += PRECISION
            z += PRECISION
        
        placed_centers = []
        for p_center in potential_centers:
            can_place = True
            for s_center in placed_centers:
                # Use squared distance for efficiency (avoids sqrt)
                dist_sq = (p_center[0] - s_center[0])**2 + \
                          (p_center[1] - s_center[1])**2 + \
                          (p_center[2] - s_center[2])**2
                if dist_sq < BALL_DIAMETER_SQ:
                    can_place = False
                    break
            if can_place:
                placed_centers.append(p_center)
                
        result = len(placed_centers)
        packer_cache[dims] = result
        return result

    # 1. Establish the baseline capacity from the original 12x12x12 box.
    TARGET_BALL_COUNT = greedy_packer(ORIGINAL_DIM, ORIGINAL_DIM, ORIGINAL_DIM)

    # 2. Search for an optimal box with integer dimensions.
    min_surface_area = ORIGINAL_SURFACE_AREA
    best_solution = None

    # Define a reasonable search range for integer dimensions (in cm).
    min_dim_search = 8
    max_dim_search = 16

    # Iterate through l, w, h where l <= w <= h to avoid permutations.
    for l in range(min_dim_search, max_dim_search + 1):
        for w in range(l, max_dim_search + 1):
            for h in range(w, max_dim_search + 1):
                
                surface_area = 2 * (l * w + l * h + w * h)
                
                # Pruning Step 1: If SA is not an improvement, skip.
                if surface_area >= min_surface_area:
                    continue
                
                # Calculate capacity. This is the most expensive step.
                num_balls = greedy_packer(float(l), float(w), float(h))
                
                # Check if this box is a better solution.
                if num_balls >= TARGET_BALL_COUNT:
                    min_surface_area = surface_area
                    best_solution = (l, w, h, surface_area)

    # 3. Output the final answer.
    if best_solution:
        l, w, h, sa = best_solution
        print(f"{l}:{w}:{h}:{sa}")
    else:
        print("0")

solve_packing_problem()