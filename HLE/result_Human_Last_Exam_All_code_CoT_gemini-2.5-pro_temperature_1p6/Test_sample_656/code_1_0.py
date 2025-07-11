import math

def solve_packing_problem():
    """
    Finds the optimal container design and ball packing to maximize energy
    subject to a surface area constraint.
    """
    # --- Constants ---
    MAX_SURFACE_AREA = 1050.0
    ENERGY_2CM_BALL = 20  # MJ
    ENERGY_1CM_BALL = 1   # MJ
    
    # --- Variables to store the best result found ---
    best_energy = 0
    best_config = {}

    # --- Search through possible box dimensions ---
    # We iterate through dimensions as integer counts of 0.5 cm units.
    # To avoid permutations, we enforce l >= w >= h (so il >= iw >= ih).
    # The loops are optimized to only check valid dimension ranges.
    
    # The maximum possible value for the smallest dimension (h) in a cube
    # with the given surface area (SA = 6*h^2)
    max_h_half_cm = int(math.sqrt(MAX_SURFACE_AREA / 6.0) / 0.5) + 1

    for ih in range(1, max_h_half_cm + 1):
        h = ih * 0.5
        
        # Max w for a given h, assuming l=w
        # SA = 2*(w^2 + 2wh) <= 1050 -> w^2 + 2wh - 525 <= 0
        # From quadratic formula: w <= -h + sqrt(h^2 + 525)
        max_w_half_cm = int((-h + math.sqrt(h*h + 525)) / 0.5) + 1
        
        for iw in range(ih, max_w_half_cm + 1):
            w = iw * 0.5

            # Max l for a given w, h
            # SA = 2*(lw + wh + hl) <= 1050 -> l(w+h) <= 525 - wh
            if (w + h) == 0: continue
            max_l_half_cm = int(((MAX_SURFACE_AREA / 2.0) - w * h) / (w + h) / 0.5) + 1
            
            for il in range(iw, max_l_half_cm + 1):
                l = il * 0.5
                
                # Final precise check on Surface Area
                surface_area = 2 * (l * w + w * h + h * l)
                if surface_area > MAX_SURFACE_AREA:
                    # Since l is increasing, further iterations in this loop will also fail
                    break

                # --- Packing Calculation ---
                # Model: Simple cubic packing, 2cm balls first, then 1cm in remainder.

                # 1. Pack 2-cm balls (diameter 4 cm)
                n2_x = math.floor(l / 4)
                n2_y = math.floor(w / 4)
                n2_z = math.floor(h / 4)
                num_2cm_balls = n2_x * n2_y * n2_z

                # 2. Pack 1-cm balls (diameter 2 cm) in the remaining space
                # The space occupied by 2cm balls forms a central box
                l_used_by_2cm = n2_x * 4
                w_used_by_2cm = n2_y * 4
                h_used_by_2cm = n2_z * 4

                # The remaining volume can be decomposed into 3 non-overlapping boxes
                rem_l = l - l_used_by_2cm
                rem_w = w - w_used_by_2cm
                rem_h = h - h_used_by_2cm
                
                # Box 1 (along length)
                n1_box1 = math.floor(rem_l / 2) * math.floor(w / 2) * math.floor(h / 2)
                # Box 2 (along width)
                n1_box2 = math.floor(l_used_by_2cm / 2) * math.floor(rem_w / 2) * math.floor(h / 2)
                # Box 3 (along height)
                n1_box3 = math.floor(l_used_by_2cm / 2) * math.floor(w_used_by_2cm / 2) * math.floor(rem_h / 2)
                
                num_1cm_balls = n1_box1 + n1_box2 + n1_box3

                # 3. Calculate total energy
                current_energy = (num_2cm_balls * ENERGY_2CM_BALL) + (num_1cm_balls * ENERGY_1CM_BALL)
                
                # 4. Check if this is the best configuration so far
                if current_energy > best_energy:
                    best_energy = current_energy
                    best_config = {
                        "type": "box",
                        "l": l, "w": w, "h": h,
                        "n1": num_1cm_balls,
                        "n2": num_2cm_balls,
                        "surface_area": surface_area
                    }

    # --- Print the final answer ---
    if best_config:
        c = f"{best_config['type']} {best_config['l']}x{best_config['w']}x{best_config['h']}"
        a = best_config['n1']
        b = best_config['n2']
        total_e = best_config['n1'] * ENERGY_1CM_BALL + best_config['n2'] * ENERGY_2CM_BALL

        print(f"Optimal design found:")
        print(f"Container: {c}")
        print(f"Number of 1-cm balls (a): {a}")
        print(f"Number of 2-cm balls (b): {b}")
        print(f"\nEnergy Calculation:")
        print(f"{ENERGY_1CM_BALL} * {a} + {ENERGY_2CM_BALL} * {b} = {total_e} MJ")
        
        final_answer = f"[{c}]{a};{b}"
        print(f"\nFinal Answer Format:")
        print(f"<<<{final_answer}>>>")

    else:
        print("No solution found within the given constraints.")
        print("<<<[0]>>>")

if __name__ == '__main__':
    solve_packing_problem()