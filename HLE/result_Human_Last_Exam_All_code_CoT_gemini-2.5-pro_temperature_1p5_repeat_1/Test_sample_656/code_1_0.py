import math

def solve_packing_problem():
    """
    Finds the optimal container design to maximize potential energy.

    This function iterates through possible box dimensions, calculates the number
    of each type of energy ball that can be packed using a grid model, and
    determines the configuration that yields the maximum total energy.
    """
    max_energy = 0
    best_config = {}
    
    # Dimensions are multiples of 0.5 cm. We iterate in steps of 0.5.
    # A generous upper bound for any dimension L is when SA = 2*L^2 = 1050 -> L ~ 22.9
    # We will search up to L=25.0 cm.
    # To handle 0.5cm increments, we can work with integers (steps) and divide by 2 later.
    # L_s = L * 2.  SA = 2(LW+WH+HL) -> SA_s = 2*(0.25*L_s*W_s + ...) -> 0.5*(L_s*W_s+...) <= 1050
    # L_s*W_s + W_s*H_s + H_s*L_s <= 2100
    
    # Loop through length, width, height in integer steps (1 step = 0.5 cm)
    # Assume L >= W >= H to avoid duplicate checks of the same box shape
    l_s_max = int(math.sqrt(2100 / 3.0)) + 2 # Upper bound for cubic shape
    
    for l_s in range(1, l_s_max):
        for w_s in range(1, l_s + 1):
            # Bound for h_s to reduce search space
            if (l_s * w_s > 2100):
                continue
            h_s_max_num = 2100 - l_s * w_s
            h_s_max_den = l_s + w_s
            if h_s_max_den == 0: continue
            h_s_max = int(h_s_max_num / h_s_max_den) + 1
            
            for h_s in range(1, min(w_s, h_s_max) + 1):
                
                # Check surface area constraint
                surface_area_s = 2 * (l_s * w_s + w_s * h_s + h_s * l_s)
                if surface_area_s > 4200: # 1050 * 4
                    continue

                # Convert steps back to cm
                L, W, H = l_s / 2.0, w_s / 2.0, h_s / 2.0

                # --- Calculate number of balls using grid packing model ---
                
                # 1. Pack large balls (diameter 4 cm)
                n2_x = math.floor(L / 4.0)
                n2_y = math.floor(W / 4.0)
                n2_z = math.floor(H / 4.0)
                n2 = n2_x * n2_y * n2_z
                
                # Dimensions of the block occupied by large balls
                Lb = n2_x * 4.0
                Wb = n2_y * 4.0
                Hb = n2_z * 4.0
                
                # 2. Pack small balls (diameter 2 cm) in the 3 remaining rectangular regions
                n1 = 0
                # Region 1: (L-Lb) x W x H
                if L - Lb >= 2:
                    n1 += math.floor((L - Lb) / 2.0) * math.floor(W / 2.0) * math.floor(H / 2.0)
                # Region 2: Lb x (W-Wb) x H
                if W - Wb >= 2:
                    n1 += math.floor(Lb / 2.0) * math.floor((W - Wb) / 2.0) * math.floor(H / 2.0)
                # Region 3: Lb x Wb x (H-Hb)
                if H - Hb >= 2:
                     n1 += math.floor(Lb / 2.0) * math.floor(Wb / 2.0) * math.floor((H - Hb) / 2.0)

                # --- Calculate total energy and update best ---
                current_energy = n2 * 20 + n1 * 1
                
                if current_energy > max_energy:
                    max_energy = current_energy
                    best_config = {
                        "type": "box",
                        "l": L, "w": W, "h": H,
                        "n1": n1,
                        "n2": n2,
                        "energy": current_energy,
                        "sa": surface_area_s / 4.0
                    }

    # Print the result
    cfg = best_config
    desc = f"{cfg['type']} {cfg['l']}x{cfg['w']}x{cfg['h']}"
    n1 = cfg['n1']
    n2 = cfg['n2']
    
    print(f"[{desc}]{int(n1)};{int(n2)}")
    print("Final equation for maximum energy:")
    print(f"{max_energy} MJ = {int(n1)} * 1 MJ + {int(n2)} * 20 MJ")

solve_packing_problem()