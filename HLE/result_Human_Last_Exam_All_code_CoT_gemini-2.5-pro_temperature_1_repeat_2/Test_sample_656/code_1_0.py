import math

def solve_pioneer_probe_packing():
    """
    Solves the Pioneer probe packing problem by searching for the optimal
    box dimensions and packing configuration.
    """
    
    # --- Configuration ---
    MAX_SURFACE_AREA = 1050.0  # cm^2
    DIM_STEP = 0.5            # cm
    
    # Energy values
    ENERGY_LARGE_BALL = 20  # MJ for 2-cm radius ball
    ENERGY_SMALL_BALL = 1   # MJ for 1-cm radius ball
    
    # Ball properties
    DIAMETER_LARGE_BALL = 4.0 # cm
    RADIUS_SMALL_BALL = 1.0   # cm
    
    # --- Search for Optimal Box Configuration ---
    
    best_config = {
        "L": 0, "W": 0, "H": 0,
        "n_large": 0, "n_small": 0,
        "energy": 0, "surface_area": 0
    }
    
    # Search ranges for dimensions (L >= W >= H to avoid redundant checks)
    # The maximum possible dimension is less than sqrt(1050/2) approx 23.
    # We start dimensions at the diameter of the large ball.
    l_range = [i * DIM_STEP for i in range(int(DIAMETER_LARGE_BALL / DIM_STEP), int(25 / DIM_STEP))]
    w_range = [i * DIM_STEP for i in range(int(DIAMETER_LARGE_BALL / DIM_STEP), int(25 / DIM_STEP))]
    h_range = [i * DIM_STEP for i in range(int(DIAMETER_LARGE_BALL / DIM_STEP), int(17 / DIM_STEP))]

    for l in l_range:
        for w in w_range:
            if w > l: continue
            for h in h_range:
                if h > w: continue
                
                # Check surface area constraint
                surface_area = 2 * (l*w + l*h + w*h)
                if surface_area > MAX_SURFACE_AREA:
                    continue
                
                # --- Calculate number of balls based on a grid packing model ---
                
                # 1. Pack large (2-cm radius) balls in a simple cubic grid
                # The number of balls that fit along a dimension is floor(Dim / Diameter)
                n_l_large = math.floor(l / DIAMETER_LARGE_BALL)
                n_w_large = math.floor(w / DIAMETER_LARGE_BALL)
                n_h_large = math.floor(h / DIAMETER_LARGE_BALL)
                num_large_balls = n_l_large * n_w_large * n_h_large
                
                if num_large_balls == 0:
                    continue

                # 2. Pack small (1-cm radius) balls in the interstitial voids
                # The voids in a simple cubic lattice are centered between the main spheres.
                # If large ball centers are at (2,6,10...), void centers are at (4,8,12...).
                # We need to find how many void centers fit within the box, leaving 1cm clearance.
                # A void center at x_c = 4+4i needs 1 <= x_c <= l-1, so 4i <= l-5. i <= (l-5)/4.
                if l >= 5 and w >= 5 and h >= 5:
                    n_l_small_voids = math.floor((l - 5) / 4) + 1
                    n_w_small_voids = math.floor((w - 5) / 4) + 1
                    n_h_small_voids = math.floor((h - 5) / 4) + 1
                    num_small_balls = n_l_small_voids * n_w_small_voids * n_h_small_voids
                else:
                    num_small_balls = 0

                # Calculate total energy
                total_energy = num_large_balls * ENERGY_LARGE_BALL + num_small_balls * ENERGY_SMALL_BALL
                
                # If this configuration is better, save it
                if total_energy > best_config["energy"]:
                    best_config = {
                        "L": l, "W": w, "H": h,
                        "n_large": num_large_balls,
                        "n_small": num_small_balls,
                        "energy": total_energy,
                        "surface_area": surface_area
                    }

    # --- Print the results ---
    b = best_config["n_large"] # b is the number of 2-cm balls
    a = best_config["n_small"] # a is the number of 1-cm balls
    L, W, H = best_config["L"], best_config["W"], best_config["H"]

    print("--- Optimal Configuration Found ---")
    print(f"Container Description (C): box {L}x{W}x{H}")
    print(f"Surface Area: {best_config['surface_area']:.2f} cm^2")
    print(f"Number of 1-cm balls (a): {a}")
    print(f"Number of 2-cm balls (b): {b}")
    print("\n--- Final Energy Calculation ---")
    print(f"Total Energy = {b} * {ENERGY_LARGE_BALL} MJ + {a} * {ENERGY_SMALL_BALL} MJ = {best_config['energy']} MJ")

    # Final answer in the required format [C]a;b
    # Note: problem statement uses 'a' for 1-cm and 'b' for 2-cm balls
    final_answer_string = f"[box {L}x{W}x{H}]{a};{b}"
    print(f"\n<<<{final_answer_string}>>>")


solve_pioneer_probe_packing()