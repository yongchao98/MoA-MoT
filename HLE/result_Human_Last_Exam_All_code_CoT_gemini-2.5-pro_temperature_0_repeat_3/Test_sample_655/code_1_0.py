import math
import numpy as np

def solve_pioneer_probe_packing():
    """
    This script finds the optimal container and ball packing to maximize energy
    for the Pioneer probe, subject to a surface area constraint.
    """
    # --- Configuration ---
    MAX_AREA = 1050.0
    STEP = 0.5
    R_LARGE = 2.0
    R_SMALL = 1.0
    E_LARGE = 10
    E_SMALL = 1

    # Use a dictionary to keep track of the best result found
    best_result = {
        "energy": -1,
        "description": "",
        "n_large": 0,
        "n_small": 0,
    }

    # --- Core Functions ---

    def distance_sq(p1, p2):
        """Calculates the squared Euclidean distance between two 3D points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def greedy_packer(is_in_bounds_func, grid_points):
        """
        Performs a greedy packing of large and small spheres.
        It prioritizes placing large balls first to maximize energy.
        """
        packed_balls = []  # List of (center_tuple, radius_float)

        # --- 1. Pack Large Balls ---
        r_large = R_LARGE
        for p in grid_points:
            if not is_in_bounds_func(p, r_large):
                continue

            can_place = True
            for packed_p, packed_r in packed_balls:
                required_dist_sq = (r_large + packed_r)**2
                if distance_sq(p, packed_p) < required_dist_sq:
                    can_place = False
                    break
            
            if can_place:
                packed_balls.append((p, r_large))

        # --- 2. Pack Small Balls in remaining space ---
        r_small = R_SMALL
        for p in grid_points:
            if not is_in_bounds_func(p, r_small):
                continue

            can_place = True
            for packed_p, packed_r in packed_balls:
                required_dist_sq = (r_small + packed_r)**2
                if distance_sq(p, packed_p) < required_dist_sq:
                    can_place = False
                    break
            
            if can_place:
                packed_balls.append((p, r_small))

        # --- 3. Calculate total energy ---
        n_large = sum(1 for _, r in packed_balls if r == R_LARGE)
        n_small = sum(1 for _, r in packed_balls if r == R_SMALL)
        energy = n_large * E_LARGE + n_small * E_SMALL
        
        return energy, n_large, n_small

    def check_and_update_best(energy, n_large, n_small, description):
        """Updates the global best_result if the current result is better."""
        nonlocal best_result
        if energy > best_result["energy"]:
            best_result = {
                "energy": energy,
                "description": description,
                "n_large": n_large,
                "n_small": n_small,
            }
            print(f"Found new best configuration: {description} with {energy} MJ")


    # --- Shape-Specific Search Functions ---

    def search_spheres():
        """Tests the largest possible sphere."""
        print("\nSearching sphere containers...")
        r = 9.0  # Largest radius with area <= 1050
        
        max_dim = int(r / STEP)
        grid_points = []
        # Generate grid points sorted from the center outwards
        for i in range(max_dim + 1):
            for j in range(max_dim + 1):
                for k in range(max_dim + 1):
                    if i**2+j**2+k**2 > max_dim**2: continue
                    # Add all 8 symmetric points
                    for sx in [-1, 1]:
                        for sy in [-1, 1]:
                            for sz in [-1, 1]:
                                if i==0 and sx==-1: continue
                                if j==0 and sy==-1: continue
                                if k==0 and sz==-1: continue
                                grid_points.append((i*STEP*sx, j*STEP*sy, k*STEP*sz))
        grid_points = sorted(list(set(grid_points)), key=lambda p: p[0]**2+p[1]**2+p[2]**2)

        def is_in_bounds(p, ball_r):
            return p[0]**2 + p[1]**2 + p[2]**2 <= (r - ball_r)**2

        energy, n_large, n_small = greedy_packer(is_in_bounds, grid_points)
        description = f"sphere r={r}"
        check_and_update_best(energy, n_large, n_small, description)

    def search_cylinders():
        """Tests a promising cylinder configuration."""
        print("\nSearching cylinder containers...")
        r, h = 7.5, 14.5 # A near-optimal configuration for volume
        
        max_r_dim = int(r / STEP)
        max_h_dim = int((h/2) / STEP)
        grid_points = []
        for k in range(max_h_dim, -1, -1):
            for i in range(max_r_dim + 1):
                for j in range(max_r_dim + 1):
                    if i**2+j**2 > max_r_dim**2: continue
                    for sx in [-1, 1]:
                        for sy in [-1, 1]:
                            for sz in [-1, 1]:
                                if i==0 and sx==-1: continue
                                if j==0 and sy==-1: continue
                                if k==0 and sz==-1: continue
                                grid_points.append((i*STEP*sx, j*STEP*sy, k*STEP*sz))
        grid_points = sorted(list(set(grid_points)), key=lambda p: (p[2]**2, p[0]**2+p[1]**2))

        def is_in_bounds(p, ball_r):
            return p[0]**2 + p[1]**2 <= (r - ball_r)**2 and abs(p[2]) <= (h/2 - ball_r)

        energy, n_large, n_small = greedy_packer(is_in_bounds, grid_points)
        description = f"cylinder r={r}, h={h}"
        check_and_update_best(energy, n_large, n_small, description)

    def search_boxes():
        """Tests a promising box configuration."""
        print("\nSearching box containers...")
        L, W, H = 13.0, 13.0, 13.0 # The largest possible cube
        
        max_l, max_w, max_h = int((L/2)/STEP), int((W/2)/STEP), int((H/2)/STEP)
        grid_points = []
        for k in range(max_h, -1, -1):
            for j in range(max_w, -1, -1):
                for i in range(max_l, -1, -1):
                    for sx in [-1, 1]:
                        for sy in [-1, 1]:
                            for sz in [-1, 1]:
                                if i==0 and sx==-1: continue
                                if j==0 and sy==-1: continue
                                if k==0 and sz==-1: continue
                                grid_points.append((i*STEP*sx, j*STEP*sy, k*STEP*sz))
        grid_points = sorted(list(set(grid_points)), key=lambda p: (p[2]**2, p[1]**2, p[0]**2))

        def is_in_bounds(p, ball_r):
            return abs(p[0]) <= (L/2 - ball_r) and \
                   abs(p[1]) <= (W/2 - ball_r) and \
                   abs(p[2]) <= (H/2 - ball_r)

        energy, n_large, n_small = greedy_packer(is_in_bounds, grid_points)
        description = f"box {L}x{W}x{H}"
        check_and_update_best(energy, n_large, n_small, description)

    # --- Main Execution ---
    search_spheres()
    search_cylinders()
    search_boxes()

    # --- Final Output ---
    b = best_result["n_large"]
    a = best_result["n_small"]
    total_energy = best_result["energy"]
    desc = best_result["description"]

    print("\n--- Final Result ---")
    print(f"The best container found is: {desc}")
    print(f"It can hold {b} large balls (2-cm radius) and {a} small balls (1-cm radius).")
    print("The final energy equation is:")
    print(f"{b} * 10 + {a} * 1 = {total_energy}")
    
    final_answer_string = f"[{desc}]{a};{b}"
    print(f"<<<{final_answer_string}>>>")

solve_pioneer_probe_packing()