import math

def solve_packing_problem():
    """
    This script finds the optimal container design (box or cylinder) to maximize
    the stored energy from two types of energy balls, subject to a surface area constraint.

    The method is an exhaustive search over possible container dimensions that are
    multiples of 0.5 cm. For each container, it calculates the maximum packable
    energy using a heuristic that prioritizes the more energy-dense large balls
    and then fills the remaining space with small balls.
    """

    # --- Configuration ---
    MAX_SA = 1050.0
    STEP = 0.5
    # Ball info:
    # Small ball: radius 1cm, energy 1MJ, diameter 2cm
    # Large ball: radius 2cm, energy 20MJ, diameter 4cm
    
    max_energy = 0
    best_config = {}

    def get_best_packing_for_box(L, W, H):
        """
        Calculates the maximum energy for a given box using a simplified
        grid-packing model that prioritizes large balls.
        """
        # A simple mixed strategy: pack large balls first in a grid,
        # then pack small balls in the three remaining disjoint rectangular slabs.
        n2_mixed = 0
        if L >= 4 and W >= 4 and H >= 4:
            n2_mixed = math.floor(L / 4) * math.floor(W / 4) * math.floor(H / 4)
        else:
            # Cannot fit any large balls, so pack only small balls.
            n1_only = math.floor(L / 2) * math.floor(W / 2) * math.floor(H / 2)
            energy = n1_only * 1
            return energy, n1_only, 0

        # Space occupied by large balls
        L_occ = 4 * math.floor(L / 4)
        W_occ = 4 * math.floor(W / 4)
        H_occ = 4 * math.floor(H / 4)

        # Calculate small balls in the 3 remaining slabs
        # Slab along L dimension
        rem_L = L - L_occ
        n1_slab_L = math.floor(rem_L / 2) * math.floor(W / 2) * math.floor(H / 2)

        # Slab along W dimension
        rem_W = W - W_occ
        n1_slab_W = math.floor(L_occ / 2) * math.floor(rem_W / 2) * math.floor(H / 2)

        # Slab along H dimension
        rem_H = H - H_occ
        n1_slab_H = math.floor(L_occ / 2) * math.floor(W_occ / 2) * math.floor(rem_H / 2)
        
        n1_mixed = n1_slab_L + n1_slab_W + n1_slab_H
        energy = n2_mixed * 20 + n1_mixed * 1

        return energy, n1_mixed, n2_mixed

    # --- Box Search ---
    # To reduce search space, assume H <= W <= L
    # For a cube, 6*L^2 <= 1050 -> L <= sqrt(175) ~= 13.2. Search up to 16 for safety.
    h_range = [i * STEP for i in range(1, int(16 / STEP) + 1)]
    
    for h in h_range:
        # Heuristic upper bound for W, assuming L=W
        w_limit = -h + math.sqrt(h*h + 525)
        w_range = [i * STEP for i in range(int(h / STEP), int(w_limit / STEP) + 2)]

        for w in w_range:
            if 2 * (w * h * 2 + w * w) > MAX_SA: continue

            # From 2(LW + WH + HL) <= 1050 -> L <= (525 - WH) / (W + H)
            if (w + h) == 0: continue
            l_limit = (525 - w * h) / (w + h)
            if l_limit < w: continue
                
            l_range = [i * STEP for i in range(int(w / STEP), int(l_limit / STEP) + 2)]
            
            for l in l_range:
                surface_area = 2 * (l * w + w * h + h * l)
                if surface_area > MAX_SA:
                    continue

                energy, n1, n2 = get_best_packing_for_box(l, w, h)
                
                if energy > max_energy:
                    max_energy = energy
                    best_config = {
                        "type": "box", "l": l, "w": w, "h": h,
                        "energy": energy, "n1": n1, "n2": n2,
                    }

    # --- Cylinder Search ---
    # 2*pi*r^2 <= 1050 -> r <= sqrt(1050/(2*pi)) ~= 12.9. Search up to 13.
    r_range = [i * STEP for i in range(1, int(13.0 / STEP) + 1)]
    
    for r in r_range:
        if 2 * math.pi * r * r > MAX_SA: continue
            
        h_limit = (MAX_SA - 2 * math.pi * r * r) / (2 * math.pi * r)
        h_range_cyl = [i * STEP for i in range(1, int(h_limit / STEP) + 2)]
        
        for h in h_range_cyl:
            surface_area = 2 * math.pi * r * r + 2 * math.pi * r * h
            if surface_area > MAX_SA: continue

            # Pack balls into the largest axis-aligned bounding box inside the cylinder
            # The side of the largest square that fits in a circle is r*sqrt(2)
            l_box = w_box = r * math.sqrt(2)
            energy, n1, n2 = get_best_packing_for_box(l_box, w_box, h)
            
            if energy > max_energy:
                max_energy = energy
                best_config = {
                    "type": "cylinder", "r": r, "h": h,
                    "energy": energy, "n1": n1, "n2": n2,
                }

    # --- Print Result ---
    if best_config:
        C_desc = ""
        if best_config["type"] == "box":
            C_desc = f"box {best_config['l']:.1f}x{best_config['w']:.1f}x{best_config['h']:.1f}"
        else: # cylinder
            C_desc = f"cylinder r={best_config['r']:.1f}, h={best_config['h']:.1f}"

        a = best_config['n1']
        b = best_config['n2']
        
        final_answer = f"[{C_desc}]{a};{b}"
        print(final_answer)
    else:
        print("[0]")

solve_packing_problem()