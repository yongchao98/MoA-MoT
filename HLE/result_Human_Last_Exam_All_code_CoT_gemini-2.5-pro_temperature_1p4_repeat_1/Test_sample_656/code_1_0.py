import math

def solve_packing_problem():
    """
    Finds the optimal container and ball packing for maximum energy.
    """
    best_energy = -1
    best_config_str = ""
    best_n1 = 0
    best_n2 = 0

    # Unit for internal calculations is 0.5 cm.
    # 2-cm radius ball (4cm diameter) -> 8x8x8 unit bounding box.
    # 1-cm radius ball (2cm diameter) -> 4x4x4 unit bounding box.

    def count_small_balls(l, w, h):
        """Counts how many 4x4x4 unit cubes fit in an l x w x h unit box."""
        if l < 4 or w < 4 or h < 4:
            return 0
        return (l // 4) * (w // 4) * (h // 4)

    def calculate_packed_energy(L, W, H):
        """
        Calculates max energy for a box of L x W x H units using a greedy strategy.
        Returns: (total_energy, num_small_balls, num_large_balls)
        """
        dims = sorted([L, W, H], reverse=True)
        L, W, H = dims[0], dims[1], dims[2]

        # Case 1: Pack only small balls (4x4x4)
        n1_only = count_small_balls(L, W, H)
        energy_only_n1 = n1_only * 1

        # Case 2: Greedy - pack large balls (8x8x8) first, then small
        if L < 8 or W < 8 or H < 8:
            # Cannot fit any large balls, so only option is small ones
            return energy_only_n1, n1_only, 0

        n2_greedy = (L // 8) * (W // 8) * (H // 8)
        energy_greedy = 20 * n2_greedy

        # Decompose remaining space into 3 non-overlapping boxes
        used_L = (L // 8) * 8
        used_W = (W // 8) * 8
        used_H = (H // 8) * 8

        n1_1 = count_small_balls(L - used_L, W, H)
        n1_2 = count_small_balls(used_L, W - used_W, H)
        n1_3 = count_small_balls(used_L, used_W, H - used_H)
        
        n1_greedy = n1_1 + n1_2 + n1_3
        energy_greedy += 1 * n1_greedy
        
        # Return the best of the two strategies
        if energy_only_n1 > energy_greedy:
            return energy_only_n1, n1_only, 0
        else:
            return energy_greedy, n1_greedy, n2_greedy

    # --- Part 1: Box Container ---
    # Surface Area (in units): LW + LH + WH <= 2100.
    # To fit at least one large ball, dimensions must be >= 8 units.
    # Assume L >= W >= H for the loops to avoid redundant shapes.
    H_unit_max = int(math.sqrt(2100 / 3))  # H <= sqrt(700) ~ 26
    for H_unit in range(8, H_unit_max + 1):
        # Loose bound for W: W*H + W*H <= 2100 -> W <= 1050/H
        W_unit_max = int(1050 / H_unit)
        for W_unit in range(H_unit, W_unit_max + 1):
            if (W_unit + H_unit) == 0: continue
            # From L(W+H) <= 2100 - WH
            L_unit_max = int((2100 - W_unit * H_unit) / (W_unit + H_unit))
            if L_unit_max < W_unit: continue
            
            # Check the box with the largest possible L for this W, H
            L_unit = L_unit_max
            
            energy, n1, n2 = calculate_packed_energy(L_unit, W_unit, H_unit)

            if energy > best_energy:
                best_energy = energy
                best_n1 = n1
                best_n2 = n2
                l, w, h = L_unit * 0.5, W_unit * 0.5, H_unit * 0.5
                best_config_str = f"box {l}x{w}x{h}"
                
    # --- Part 2: Cylinder Container ---
    # Surface Area (in cm): r^2 + rh <= 1050 / (2*pi) ~ 167.1
    # To fit a large ball, inscribed box needs side >= 4cm, so sqrt(2)*r >= 4 -> r >= 2.83cm.
    # And height must be >= 4cm.
    for r_cm in [i * 0.5 for i in range(int(3.0*2), int(12.9*2) + 1)]: # r from 3.0 to 12.5
        h_cm_max = (167.1 - r_cm**2) / r_cm
        if h_cm_max < 4.0: continue
        for h_cm in [i * 0.5 for i in range(int(4.0*2), int(h_cm_max*2) + 1)]: # h from 4.0
            # Inscribe a square-base box inside the cylinder
            R_unit, H_unit = int(r_cm * 2), int(h_cm * 2)
            box_side_unit = int(math.sqrt(2) * R_unit)
            
            energy, n1, n2 = calculate_packed_energy(box_side_unit, box_side_unit, H_unit)
            
            if energy > best_energy:
                best_energy = energy
                best_n1 = n1
                best_n2 = n2
                best_config_str = f"cylinder r={r_cm}, h={h_cm}"

    print(f"[{best_config_str}]{best_n1};{best_n2}")

solve_packing_problem()