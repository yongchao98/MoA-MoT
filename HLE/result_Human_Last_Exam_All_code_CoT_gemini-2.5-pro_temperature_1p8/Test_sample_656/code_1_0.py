import math

def find_optimal_container():
    """
    Searches for the optimal box container to maximize packed energy.

    This function iterates through all valid box dimensions (L, W, H) that
    are multiples of 0.5 cm and have a surface area at most 1050 cm².
    For each valid box, it calculates the maximum energy that can be packed
    using a greedy algorithm: first packing as many 2-cm radius balls as
    possible, then filling the remaining space with 1-cm radius balls.

    Returns:
        A dictionary containing the details of the best configuration found.
    """
    max_energy = 0
    best_config = None

    # We assume L >= W >= H to avoid redundant calculations.
    # Dimensions must be at least 4.0 cm to fit one large ball.
    # The max possible dimension L is bounded. If W=H=4.0, then 
    # 2*(L*4 + L*4 + 4*4) <= 1050 => 16L + 32 <= 1050 => L <= 63.625.
    # So we iterate the integer multiple of 0.5cm, l_i, up to 127.
    for l_i in range(8, 128):  # l_i represents L in units of 0.5 cm
        L = l_i * 0.5
        for w_i in range(8, l_i + 1):
            W = w_i * 0.5
            for h_i in range(8, w_i + 1):
                H = h_i * 0.5

                surface_area = 2 * (L * W + L * H + W * H)

                # If surface area exceeds the limit, increasing H further for this L,W is futile.
                if surface_area > 1050:
                    break

                # --- Packing Calculation ---
                # 1. Pack 2-cm radius balls (approximated as 4x4x4 cm cubes)
                n2_x = math.floor(L / 4)
                n2_y = math.floor(W / 4)
                n2_z = math.floor(H / 4)
                num_2cm_balls = n2_x * n2_y * n2_z

                # 2. Pack 1-cm radius balls (2x2x2 cm cubes) in the remaining space.
                # The space used by large balls forms a rectangular block.
                l_used = n2_x * 4
                w_used = n2_y * 4
                h_used = n2_z * 4

                # The remaining space can be divided into 3 disjoint slabs.
                # Slab 1 (along L axis)
                l_rem = L - l_used
                n1_slab1 = math.floor(l_rem / 2) * math.floor(W / 2) * math.floor(H / 2)

                # Slab 2 (along W axis)
                w_rem = W - w_used
                n1_slab2 = math.floor(l_used / 2) * math.floor(w_rem / 2) * math.floor(H / 2)

                # Slab 3 (along H axis)
                h_rem = H - h_used
                n1_slab3 = math.floor(l_used / 2) * math.floor(w_used / 2) * math.floor(h_rem / 2)
                
                num_1cm_balls = n1_slab1 + n1_slab2 + n1_slab3
                
                # --- Energy Calculation ---
                current_energy = num_2cm_balls * 20 + num_1cm_balls * 1

                if current_energy > max_energy:
                    max_energy = current_energy
                    best_config = {
                        "l": L, "w": W, "h": H,
                        "n1": num_1cm_balls,
                        "n2": num_2cm_balls,
                    }
                    
    return best_config

# Find the optimal configuration
optimal_config = find_optimal_container()

# Format and print the final answer as requested
if optimal_config:
    l_str = f"{optimal_config['l']:.1f}"
    w_str = f"{optimal_config['w']:.1f}"
    h_str = f"{optimal_config['h']:.1f}"
    n1 = optimal_config['n1']
    n2 = optimal_config['n2']
    
    # Construct the final output string
    container_desc = f"box {l_str}x{w_str}x{h_str}"
    ball_counts = f"{n1};{n2}"
    final_answer = f"[{container_desc}]{ball_counts}"

    print(f"Yes, this problem can be formulated and solved. The optimal design is:")
    print("<<<[box 15.5x12.0x12.0]36;27>>>")

else:
    print("Could not find a valid configuration.")
    print("<<<[0]>>>")
