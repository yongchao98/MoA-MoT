import math

def solve_pioneer_probe_packing():
    """
    Finds the optimal box dimensions and ball quantities to maximize the packed potential energy.
    The method iterates through all valid box dimensions, calculates the number of packing "cells"
    for each, and then determines the best combination of energy balls based on a greedy strategy.
    """
    max_energy = 0
    best_config = {}

    # Dimensions l, w, h are in 0.5 cm units to ensure they are multiples of 0.5.
    # The surface area constraint is 2*(L*W + L*H + W*H) <= 1050.
    # In our integer units, this is 2 * (0.25*l*w + 0.25*l*h + 0.25*w*h) <= 1050
    # which simplifies to: l*w + l*h + w*h <= 2100.

    # To avoid duplicate permutations of dimensions, we enforce l <= w <= h.
    # This gives us limits for the loops. From 3*l^2 <= 2100, we get l <= sqrt(700) ~= 26.
    l_limit = int((2100 / 3)**0.5)
    for l in range(1, l_limit + 1):
        
        # From l*w + l*h + w*h >= l*w + l*w + w*w (since h>=w), we find limits for w
        w_limit = int(((2100 - l*l)/(2*l))) if (2*l) > 0 else l
        for w in range(l, w_limit + 1):
            
            # From l*w + h*(l+w) <= 2100, we find the limit for h.
            if (l+w) == 0: continue
            h_limit = int((2100 - l*w)/(l+w))
            for h in range(w, h_limit + 1):
                L, W, H = l / 2.0, w / 2.0, h / 2.0

                # A ball needs space at least its diameter (2cm for smallest ball).
                if L < 2.0 or W < 2.0 or H < 2.0:
                    continue
                
                # Calculate the number of 2x2x2 cm cells that fit in the box.
                num_cells = math.floor(L / 2) * math.floor(W / 2) * math.floor(H / 2)

                # Determine the number of each ball type.
                # n_A for 2-cm radius balls (Type A), n_B for 1-cm radius balls (Type B).
                n_A = num_cells // 8
                n_B = num_cells % 8

                current_energy = 10 * n_A + n_B
                
                if current_energy > max_energy:
                    max_energy = current_energy
                    best_config = {
                        "l": L, "w": W, "h": H,
                        "n_A": n_A, # 2-cm balls
                        "n_B": n_B  # 1-cm balls
                    }

    # Format the output as specified in the problem description.
    # [C]a;b where a is number of 1-cm balls, b is number of 2-cm balls.
    l, w, h = best_config['l'], best_config['w'], best_config['h']
    
    # Use integer dimensions if they are whole numbers for cleaner output
    l_str = str(int(l)) if l == int(l) else str(l)
    w_str = str(int(w)) if w == int(w) else str(w)
    h_str = str(int(h)) if h == int(h) else str(h)
    
    container_desc = f"box {l_str}x{w_str}x{h_str}"
    
    print(f"[{container_desc}]{best_config['n_B']};{best_config['n_A']}")

solve_pioneer_probe_packing()
<<<box 10x14x16;0;35>>>