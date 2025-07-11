import math

def solve_packing_problem():
    """
    Finds the optimal container design to maximize potential energy.

    This function searches over possible box and cylinder dimensions to find the
    configuration that yields the maximum energy from packed spherical balls,
    subject to a surface area constraint. It assumes a simplified cubic
    packing model for tractability.
    """
    max_energy = 0.0
    best_config = {}
    MAX_SA = 1050.0
    
    # --- 1. Box Container Search ---
    # Search range for dimensions (multiples of 0.5). A dimension larger than ~25 cm
    # is unlikely given the SA constraint.
    dim_range = [i * 0.5 for i in range(1, 51)] 
    
    for l in dim_range:
        for w in dim_range:
            # To reduce redundant calculations (e.g., 8x16x16 vs 16x8x16)
            if w < l:
                continue
            for h in dim_range:
                if h < w:
                    continue
                
                sa = 2 * (l*w + l*h + w*h)
                if sa <= MAX_SA:
                    # Calculate number of 2-cm balls (in 4x4x4 cm cubes)
                    n2_l = math.floor(l / 4)
                    n2_w = math.floor(w / 4)
                    n2_h = math.floor(h / 4)
                    n2 = n2_l * n2_w * n2_h
                    
                    # Calculate remaining space for 1-cm balls (in 2x2x2 cm cubes)
                    l_used_by_n2 = n2_l * 4
                    w_used_by_n2 = n2_w * 4
                    h_used_by_n2 = n2_h * 4
                    
                    # Leftover space is comprised of 3 slabs
                    l_rem = l - l_used_by_n2
                    w_rem = w - w_used_by_n2
                    h_rem = h - h_used_by_n2
                    
                    # Slab 1
                    n1_s1 = math.floor(l_rem / 2) * math.floor(w / 2) * math.floor(h_used_by_n2 / 2)
                    # Slab 2
                    n1_s2 = math.floor(l / 2) * math.floor(w_rem / 2) * math.floor(h_used_by_n2 / 2)
                    # Slab 3
                    n1_s3 = math.floor(l / 2) * math.floor(w / 2) * math.floor(h_rem / 2)
                    n1 = n1_s1 + n1_s2 + n1_s3

                    energy = 20 * n2 + 1 * n1
                    if energy > max_energy:
                        max_energy = energy
                        best_config = {
                            "type": "box",
                            "dims": sorted([l, w, h]),
                            "n1": n1,
                            "n2": n2,
                            "energy": energy,
                            "sa": sa
                        }

    # --- 2. Cylinder Container Search ---
    r_range = [i * 0.5 for i in range(1, 31)] # r up to 15
    h_range = [i * 0.5 for i in range(1, 61)] # h up to 30

    for r in r_range:
        for h in h_range:
            sa = 2 * math.pi * r * (r + h)
            if sa <= MAX_SA:
                # Pack inside the largest inscribed square prism
                s = r * math.sqrt(2) # Side of inscribed square
                
                # N2 balls
                n2_s = math.floor(s / 4)
                n2_h = math.floor(h / 4)
                n2 = n2_s * n2_s * n2_h

                # N1 balls in leftover height slab
                s_used_by_n2 = n2_s * 4
                h_rem = h - (n2_h * 4)
                
                n1 = math.floor(s_used_by_n2 / 2) * math.floor(s_used_by_n2 / 2) * math.floor(h_rem / 2)

                energy = 20 * n2 + 1 * n1
                if energy > max_energy:
                    max_energy = energy
                    best_config = {
                        "type": "cylinder",
                        "dims": [r, h],
                        "n1": n1,
                        "n2": n2,
                        "energy": energy,
                        "sa": sa
                    }
    
    # --- 3. Format and Print Result ---
    if best_config["type"] == "box":
        dims = best_config["dims"]
        desc = f'box {dims[0]}x{dims[1]}x{dims[2]}'
    else: # cylinder
        r, h = best_config["dims"]
        desc = f'cylinder r={r}, h={h}'
        
    n1 = best_config["n1"]
    n2 = best_config["n2"]

    final_answer = f'[{desc}]{int(n1)};{int(n2)}'
    
    print(final_answer)


solve_packing_problem()
