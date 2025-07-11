import math

def solve_packing_problem():
    """
    Finds the optimal box container and ball packing to maximize energy
    subject to a surface area constraint using a heuristic search.
    """
    max_surface_area = 1050.0
    precision = 0.5
    
    # Ball properties
    # Type 1: 1cm radius (2cm diameter), 1 MJ
    # Type 2: 2cm radius (4cm diameter), 20 MJ
    d1 = 2.0
    e1 = 1
    d2 = 4.0
    e2 = 20

    max_energy = 0
    best_config = {}

    # Define a reasonable search range for dimensions.
    # A cube with SA=1050 would have a side of ~13.2 cm.
    # We search a bit beyond that to allow for non-cubic shapes.
    min_dim_steps = int(d2 / precision)
    max_dim_steps = int(25.0 / precision)

    # Create a list of possible dimension values
    dim_values = [i * precision for i in range(min_dim_steps, max_dim_steps + 1)]

    # Iterate through all possible box dimensions (L x W x H)
    # To avoid redundant permutations (e.g., 10x12x15 vs 12x10x15),
    # we enforce L >= W >= H.
    for l_val in reversed(dim_values):
        for w_val in reversed(dim_values):
            if w_val > l_val:
                continue
            for h_val in reversed(dim_values):
                if h_val > w_val:
                    continue

                # The smallest dimension must be able to fit the largest ball
                if h_val < d2:
                    continue

                # Check surface area constraint
                surface_area = 2 * (l_val * w_val + l_val * h_val + w_val * h_val)
                if surface_area > max_surface_area:
                    continue

                # --- Start of Greedy Packing Heuristic ---

                # 1. Pack with large balls (diameter d2) first
                nx = math.floor(l_val / d2)
                ny = math.floor(w_val / d2)
                nz = math.floor(h_val / d2)
                num_b2 = nx * ny * nz

                # 2. Decompose remaining space and pack with small balls (diameter d1)
                rem_l = l_val - nx * d2
                rem_w = w_val - ny * d2
                rem_h = h_val - nz * d2

                # Calculate number of small balls that can fit in each of the 7 remaining sub-volumes
                # (The 8th sub-volume is the core, already filled with large balls)
                
                # Dimensions of the core block (filled with large balls), divided by d1
                core_l_div_d1 = math.floor(nx * d2 / d1)
                core_w_div_d1 = math.floor(ny * d2 / d1)
                core_h_div_d1 = math.floor(nz * d2 / d1)

                # Dimensions of the remainder slabs, divided by d1
                rem_l_div_d1 = math.floor(rem_l / d1)
                rem_w_div_d1 = math.floor(rem_w / d1)
                rem_h_div_d1 = math.floor(rem_h / d1)

                num_b1 = 0
                # Slab along L dimension
                num_b1 += rem_l_div_d1 * core_w_div_d1 * core_h_div_d1
                # Slab along W dimension
                num_b1 += core_l_div_d1 * rem_w_div_d1 * core_h_div_d1
                # Slab along H dimension
                num_b1 += core_l_div_d1 * core_w_div_d1 * rem_h_div_d1
                # Edge between L and W slabs
                num_b1 += rem_l_div_d1 * rem_w_div_d1 * core_h_div_d1
                # Edge between L and H slabs
                num_b1 += rem_l_div_d1 * core_w_div_d1 * rem_h_div_d1
                # Edge between W and H slabs
                num_b1 += core_l_div_d1 * rem_w_div_d1 * rem_h_div_d1
                # Corner block
                num_b1 += rem_l_div_d1 * rem_w_div_d1 * rem_h_div_d1
                
                # --- End of Heuristic ---

                current_energy = num_b1 * e1 + num_b2 * e2

                if current_energy > max_energy:
                    max_energy = current_energy
                    best_config = {
                        "type": "box",
                        "dims": (l_val, w_val, h_val),
                        "n1": num_b1,
                        "n2": num_b2,
                        "energy": current_energy,
                        "sa": surface_area
                    }

    # Print the results
    if best_config:
        dims = best_config["dims"]
        n1 = best_config["n1"]
        n2 = best_config["n2"]
        energy = best_config["energy"]
        
        desc = f"box {dims[0]}x{dims[1]}x{dims[2]}"
        final_answer_str = f"[{desc}]{n1};{n2}"
        
        print(f"Optimal container found: {desc}")
        print(f"Surface Area: {best_config['sa']:.2f} cm2")
        print(f"Packed with {n1} 1-cm balls and {n2} 2-cm balls.")
        print(f"Maximum energy calculation:")
        print(f"{n1} * {e1} MJ + {n2} * {e2} MJ = {energy} MJ")
        print("\nFinal Answer Format:")
        print(f"<<<{final_answer_str}>>>")
    else:
        print("No valid configuration found.")
        print("<<<[0]>>>")

solve_packing_problem()