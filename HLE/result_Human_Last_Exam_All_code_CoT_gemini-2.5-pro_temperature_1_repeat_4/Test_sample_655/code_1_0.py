import math
import numpy as np

def solve_packing_problem():
    """
    This script finds the optimal container to maximize stored energy.
    It works as follows:
    1.  Define three container shapes: Sphere, Cylinder, and Box.
    2.  For each shape, iterate through all possible dimensions that satisfy the surface area constraint (<= 1050 cm^2) and precision (multiples of 0.5 cm).
    3.  For each valid container, use a greedy packing algorithm to fill it with energy balls.
        a. The algorithm first places the more valuable 2-cm radius balls. It tries all possible grid locations, sorted from bottom to top, and places a ball if it doesn't overlap with others.
        b. Then, it fills the remaining gaps with 1-cm radius balls using the same method.
    4.  The total energy for the packed container is calculated (10 MJ per 2-cm ball, 1 MJ per 1-cm ball).
    5.  The script keeps track of the best configuration found (container shape, dimensions, and ball counts) that yields the highest energy.
    6.  Finally, it prints the details of the best configuration found.
    Note: This computation can be intensive and may take a few minutes to run.
    """

    # --- Configuration ---
    PRECISION = 0.5
    MAX_SA = 1050.0
    RADIUS1, ENERGY1 = 1.0, 1
    RADIUS2, ENERGY2 = 2.0, 10

    # --- Global state to store the best result ---
    best_solution = {
        "energy": -1,
        "shape": "",
        "dims_str": "",
        "n1": 0,
        "n2": 0,
    }

    # --- Helper functions ---
    def round_to_precision(value):
        return round(value / PRECISION) * PRECISION

    def distance_sq(p1, p2):
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    # --- Geometry-specific functions to get valid center locations ---
    def get_candidate_centers(shape, dims, radius):
        centers = []
        if shape == "box":
            L, W, H = dims
            for x in np.arange(radius, L - radius + 1e-9, PRECISION):
                for y in np.arange(radius, W - radius + 1e-9, PRECISION):
                    for z in np.arange(radius, H - radius + 1e-9, PRECISION):
                        centers.append((x, y, z))
        elif shape == "sphere":
            R = dims
            center_point = (R, R, R)
            max_dist_sq = (R - radius)**2
            for x in np.arange(radius, 2 * R - radius + 1e-9, PRECISION):
                for y in np.arange(radius, 2 * R - radius + 1e-9, PRECISION):
                    for z in np.arange(radius, 2 * R - radius + 1e-9, PRECISION):
                        if distance_sq((x, y, z), center_point) <= max_dist_sq + 1e-9:
                            centers.append((x, y, z))
        elif shape == "cylinder":
            r, h = dims
            center_axis_point = (r, r)
            max_radial_dist_sq = (r - radius)**2
            for x in np.arange(radius, 2 * r - radius + 1e-9, PRECISION):
                for y in np.arange(radius, 2 * r - radius + 1e-9, PRECISION):
                    if (x - center_axis_point[0])**2 + (y - center_axis_point[1])**2 <= max_radial_dist_sq + 1e-9:
                        for z in np.arange(radius, h - radius + 1e-9, PRECISION):
                            centers.append((x, y, z))
        return centers

    # --- Main packing algorithm ---
    def pack_container(shape, dims):
        placed_balls = []

        # Phase 1: Pack 2-cm balls
        candidate_centers_2 = get_candidate_centers(shape, dims, RADIUS2)
        candidate_centers_2.sort(key=lambda p: (p[2], p[1], p[0])) # Pack layer by layer
        min_dist_sq_22 = (RADIUS2 + RADIUS2)**2
        for center in candidate_centers_2:
            if all(distance_sq(center, p_center) >= min_dist_sq_22 - 1e-9 for p_center, _ in placed_balls):
                placed_balls.append((center, RADIUS2))
        n2 = len(placed_balls)

        # Phase 2: Pack 1-cm balls
        candidate_centers_1 = get_candidate_centers(shape, dims, RADIUS1)
        candidate_centers_1.sort(key=lambda p: (p[2], p[1], p[0]))
        min_dist_sq_11 = (RADIUS1 + RADIUS1)**2
        min_dist_sq_12 = (RADIUS1 + RADIUS2)**2
        n1_added = 0
        for center in candidate_centers_1:
            can_place = True
            for p_center, p_radius in placed_balls:
                required_dist_sq = min_dist_sq_12 if p_radius == RADIUS2 else min_dist_sq_11
                if distance_sq(center, p_center) < required_dist_sq - 1e-9:
                    can_place = False
                    break
            if can_place:
                placed_balls.append((center, RADIUS1))
                n1_added += 1
        n1 = n1_added
        
        energy = n1 * ENERGY1 + n2 * ENERGY2
        return energy, n1, n2

    # --- Search loop for all containers ---
    def check_and_update(energy, n1, n2, shape, dims_str):
        nonlocal best_solution
        if energy > best_solution["energy"]:
            best_solution.update({"energy": energy, "shape": shape, "dims_str": dims_str, "n1": n1, "n2": n2})

    # 1. Sphere
    max_R = math.sqrt(MAX_SA / (4 * math.pi))
    for R_float in np.arange(PRECISION, max_R + 1e-9, PRECISION):
        R = round_to_precision(R_float)
        if 4 * math.pi * R**2 > MAX_SA: continue
        energy, n1, n2 = pack_container("sphere", R)
        check_and_update(energy, n1, n2, "sphere", f"r={R}")

    # 2. Cylinder
    max_r_cyl = math.sqrt(MAX_SA / (2 * math.pi))
    for r_float in np.arange(PRECISION, max_r_cyl + 1e-9, PRECISION):
        r = round_to_precision(r_float)
        if 2 * math.pi * r**2 > MAX_SA: continue
        max_h = (MAX_SA - 2 * math.pi * r**2) / (2 * math.pi * r)
        h = round_to_precision(max_h)
        if h < PRECISION: continue
        if 2 * math.pi * r**2 + 2 * math.pi * r * h > MAX_SA: continue
        energy, n1, n2 = pack_container("cylinder", (r, h))
        check_and_update(energy, n1, n2, "cylinder", f"r={r}, h={h}")

    # 3. Box
    max_H_box = math.sqrt(MAX_SA / 6)
    for H_float in np.arange(PRECISION, max_H_box + 1e-9, PRECISION):
        H = round_to_precision(H_float)
        max_W_for_H = -H + math.sqrt(H**2 + MAX_SA / 2)
        for W_float in np.arange(H, max_W_for_H + 1e-9, PRECISION):
            W = round_to_precision(W_float)
            if 2 * H * W >= MAX_SA: continue
            max_L_num = MAX_SA - 2 * H * W
            max_L_den = 2 * (H + W)
            if max_L_den <= 0: continue
            max_L = max_L_num / max_L_den
            if max_L < W: continue
            L = round_to_precision(max_L)
            if L < W: continue
            if 2 * (L*W + L*H + W*H) > MAX_SA: continue
            energy, n1, n2 = pack_container("box", (L, W, H))
            check_and_update(energy, n1, n2, "box", f"{L}x{W}x{H}")

    # --- Print final result ---
    bs = best_solution
    container_str = f"{bs['shape']} {bs['dims_str']}"
    print(f"Optimal container found: {container_str}")
    print(f"Number of 1-cm balls (a): {bs['n1']}")
    print(f"Number of 2-cm balls (b): {bs['n2']}")
    print(f"Maximum potential energy: {bs['energy']} MJ")
    print(f"Final Answer: [{container_str}]{bs['n1']};{bs['n2']}")

solve_packing_problem()