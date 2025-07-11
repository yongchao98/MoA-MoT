import math

def get_n_squares_in_circle(R, s):
    """
    Calculates how many non-overlapping squares of side 's' can be packed
    into a circle of radius 'R', assuming the squares' centers are on a grid with spacing 's'.
    """
    count = 0
    hs = s / 2.0  # half-side of the square
    max_idx = math.floor(R / s)
    
    # Iterate through grid positions for the center of the square
    for i in range(-max_idx, max_idx + 1):
        for j in range(-max_idx, max_idx + 1):
            cx = i * s
            cy = j * s
            # Check if the corner of the square furthest from the circle's center is inside
            if (abs(cx) + hs)**2 + (abs(cy) + hs)**2 <= R**2:
                count += 1
    return count

def solve():
    """
    Finds the optimal container design by iterating through possible dimensions.
    """
    best_energy = -1
    best_config = {}
    MAX_SA = 1050.0

    # --- Part 1: Search for the best Box container ---
    # We use integers for looping (1 unit = 0.5 cm) to avoid float precision issues.
    # l, w, h are dimensions multiplied by 2.
    # SA constraint: 2(LW+LH+WH) <= 1050  =>  lw+lh+wh <= 2100
    # Assume l >= w >= h to avoid redundant checks.
    h_int_max = int((2100 / 3)**0.5)
    for h_int in range(1, h_int_max + 2):
        H = h_int / 2.0
        w_int_max = int(-h_int + (h_int**2 + 2100)**0.5)
        for w_int in range(h_int, w_int_max + 2):
            W = w_int / 2.0
            if (w_int + h_int) == 0: continue
            l_int_max = int((2100 - w_int * h_int) / (w_int + h_int))
            for l_int in range(w_int, l_int_max + 2):
                L = l_int / 2.0
                
                if 2 * (L*W + L*H + W*H) > MAX_SA:
                    continue

                n2 = math.floor(L/4) * math.floor(W/4) * math.floor(H/4)
                total_n1_slots = math.floor(L/2) * math.floor(W/2) * math.floor(H/2)
                n1 = total_n1_slots - n2 * 8
                energy = 20 * n2 + 1 * n1
                
                if energy > best_energy:
                    best_energy = energy
                    best_config = {"type": "box", "L": L, "W": W, "H": H, "n1": n1, "n2": n2}

    # --- Part 2: Search for the best Cylinder container ---
    # SA constraint: 2*pi*R^2 + 2*pi*R*H <= 1050
    SA_CONST = MAX_SA / (2 * math.pi) # R^2 + RH <= SA_CONST
    r_int_max = int((SA_CONST**0.5) * 2)
    for r_int in range(1, r_int_max + 2):
        R = r_int / 2.0
        if R == 0: continue
        h_max = (SA_CONST - R**2) / R
        if h_max < 0: continue
        h_int_max = int(h_max * 2)
        
        for h_int in range(1, h_int_max + 2):
            H = h_int / 2.0
            
            if (2 * math.pi * R**2 + 2 * math.pi * R * H) > MAX_SA:
                continue
                
            n2_per_layer = get_n_squares_in_circle(R, 4.0)
            n1_slots_per_layer = get_n_squares_in_circle(R, 2.0)
            
            n2_layers = math.floor(H / 4.0)
            n1_layers = math.floor(H / 2.0)
            
            n2 = n2_per_layer * n2_layers
            total_n1_slots = n1_slots_per_layer * n1_layers
            n1 = total_n1_slots - n2 * 8
            energy = 20 * n2 + 1 * n1
            
            if energy > best_energy:
                best_energy = energy
                best_config = {"type": "cylinder", "R": R, "H": H, "n1": n1, "n2": n2}
    
    # --- Part 3: Output the result ---
    n1 = best_config["n1"]
    n2 = best_config["n2"]
    final_energy = 20 * n2 + 1 * n1

    if best_config["type"] == "box":
        desc = "box {:g}x{:g}x{:g}".format(best_config['L'], best_config['W'], best_config['H'])
    else:
        desc = "cylinder r={:g}, h={:g}".format(best_config['R'], best_config['H'])

    print(f"The optimal container configuration is a {best_config['type']}.")
    print(f"Description: {desc}")
    print(f"Number of 1-cm balls (a): {n1}")
    print(f"Number of 2-cm balls (b): {n2}")
    print(f"\nThis configuration yields a total energy of {final_energy} MJ, based on the calculation:")
    print(f"{final_energy} = 20 * {n2} + 1 * {n1}")
    
    return f"[{desc}]{n1};{n2}"

# Run the solver and print the final answer in the required format
final_answer = solve()
print(f"\n<<< {final_answer} >>>")