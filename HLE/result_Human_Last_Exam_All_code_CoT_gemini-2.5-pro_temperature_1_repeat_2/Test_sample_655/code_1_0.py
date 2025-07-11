import math
import numpy as np

def solve_packing_problem():
    """
    This script finds the optimal container and ball configuration to maximize stored energy.
    It simulates packing for three high-potential container shapes and outputs the best result.
    """

    # --- Configuration ---
    MAX_AREA = 1050.0
    STEP = 0.5
    R_LARGE = 2.0
    R_SMALL = 1.0
    E_LARGE = 10
    E_SMALL = 1

    # --- Helper Functions ---
    def check_overlap(center1, r1, ball_list):
        """Checks if a new ball at center1 with radius r1 overlaps with any ball in ball_list."""
        for ball in ball_list:
            # Using squared distances to avoid costly sqrt
            dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(center1, ball['center']))
            if dist_sq < (r1 + ball['r'])**2:
                return True
        return False

    def is_in_container(center, radius, container_type, dims):
        """Checks if a ball is fully inside the container boundaries."""
        cx, cy, cz = center
        if container_type == "sphere":
            R = dims['R']
            # Container's bounding box is from (0,0,0) to (2R,2R,2R). Sphere center is at (R,R,R).
            return math.sqrt((cx - R)**2 + (cy - R)**2 + (cz - R)**2) + radius <= R
        elif container_type == "box":
            l, w, h = dims['l'], dims['w'], dims['h']
            # Container is a box from (0,0,0) to (l,w,h).
            return (radius <= cx <= l - radius and
                    radius <= cy <= w - radius and
                    radius <= cz <= h - radius)
        elif container_type == "cylinder":
            r, h = dims['r']
            # Container's bounding box from (0,0,0) to (2r,2r,h). Cylinder base center is at (r,r,0).
            return (math.sqrt((cx - r)**2 + (cy - r)**2) + radius <= r and
                    radius <= cz <= h - radius)
        return False

    def pack_container(container_type, dims):
        """
        Performs a greedy packing simulation for a given container.
        It first packs large balls, then small balls.
        """
        placed_balls = []
        
        # Define bounding box for the iteration grid
        if container_type == "box":
            l, w, h = dims['l'], dims['w'], dims['h']
            bb_x, bb_y, bb_z = l, w, h
        elif container_type == "sphere":
            R = dims['R']
            bb_x, bb_y, bb_z = 2 * R, 2 * R, 2 * R
        elif container_type == "cylinder":
            r, h = dims['r'], dims['h']
            bb_x, bb_y, bb_z = 2 * r, 2 * r, h

        # --- Pack large balls (radius = 2.0 cm) ---
        radius = R_LARGE
        # Iterate through possible center locations on the 0.5cm grid
        # We iterate from top-to-bottom (like gravity) for potentially better packing
        z_coords = np.arange(bb_z - radius, radius - STEP, -STEP)
        y_coords = np.arange(radius, bb_y - radius + STEP, STEP)
        x_coords = np.arange(radius, bb_x - radius + STEP, STEP)
        
        for cz in z_coords:
            for cy in y_coords:
                for cx in x_coords:
                    center = (cx, cy, cz)
                    if is_in_container(center, radius, container_type, dims):
                        if not check_overlap(center, radius, placed_balls):
                            placed_balls.append({'center': center, 'r': radius})

        # --- Pack small balls (radius = 1.0 cm) in the remaining space ---
        num_large_balls = len(placed_balls)
        radius = R_SMALL
        z_coords = np.arange(bb_z - radius, radius - STEP, -STEP)
        y_coords = np.arange(radius, bb_y - radius + STEP, STEP)
        x_coords = np.arange(radius, bb_x - radius + STEP, STEP)

        for cz in z_coords:
            for cy in y_coords:
                for cx in x_coords:
                    center = (cx, cy, cz)
                    if is_in_container(center, radius, container_type, dims):
                        if not check_overlap(center, radius, placed_balls):
                            placed_balls.append({'center': center, 'r': radius})
        
        num_small_balls = len(placed_balls) - num_large_balls
        energy = num_large_balls * E_LARGE + num_small_balls * E_SMALL
        
        return num_small_balls, num_large_balls, energy

    # --- Main Execution ---
    
    # Define the three most promising container candidates
    candidates = [
        {"shape": "sphere", "dims": {"R": 9.0}, "desc": "sphere r=9.0"},
        {"shape": "cylinder", "dims": {"r": 7.5, "h": 14.5}, "desc": "cylinder r=7.5, h=14.5"},
        {"shape": "box", "dims": {"l": 13.0, "w": 13.0, "h": 13.0}, "desc": "box 13.0x13.0x13.0"}
    ]

    best_result = {"energy": -1}

    print("Starting packing simulation for candidate containers...")
    for cand in candidates:
        print(f"\nTesting container: {cand['desc']}")
        n_small, n_large, energy = pack_container(cand['shape'], cand['dims'])
        print(f"Result: Found space for {n_large} large balls and {n_small} small balls.")
        print(f"Total Energy Calculation: {n_large} * {E_LARGE} MJ + {n_small} * {E_SMALL} MJ = {energy} MJ")

        if energy > best_result["energy"]:
            best_result = {
                "energy": energy,
                "shape_desc": cand['desc'],
                "n_small": n_small,
                "n_large": n_large
            }

    print("\n-------------------------------------------")
    print("           OPTIMAL CONFIGURATION           ")
    print("-------------------------------------------")
    print(f"Container: {best_result['shape_desc']}")
    print(f"Number of 1-cm balls (1 MJ): {best_result['n_small']}")
    print(f"Number of 2-cm balls (10 MJ): {best_result['n_large']}")
    print(f"Maximum Energy: {best_result['energy']} MJ")
    
    # Final answer in the required format [C]a;b
    final_answer = f"[{best_result['shape_desc}}]{best_result['n_small']};{best_result['n_large']}"
    print(f"\nFinal Answer String: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_packing_problem()