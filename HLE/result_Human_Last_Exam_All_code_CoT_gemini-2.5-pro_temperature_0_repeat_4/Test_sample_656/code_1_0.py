import math

def solve_pioneer_probe_packing():
    """
    Solves the energy ball packing problem for the Pioneer probe.

    This function searches for the optimal container (box or cylinder) with a
    surface area at most 1050 cm^2 to maximize the total potential energy
    from packed energy balls.

    The search iterates through possible dimensions (multiples of 0.5 cm) for
    both container types and calculates the maximum packable energy for each,
    assuming a grid-based packing strategy.
    """
    max_energy = -1
    best_config = None
    max_surface_area = 1050.0

    def calculate_energy(L, W, H):
        """
        Calculates the maximum energy packable in a box of given dimensions.
        
        This assumes a grid packing strategy, prioritizing the more energy-dense
        2-cm radius balls.
        """
        if L <= 0 or W <= 0 or H <= 0:
            return 0, 0, 0

        # Number of 2-cm radius balls (4x4x4 cm bounding box)
        n2 = math.floor(L / 4) * math.floor(W / 4) * math.floor(H / 4)
        
        # Total number of 1-cm radius ball slots (2x2x2 cm bounding box)
        total_small_slots = math.floor(L / 2) * math.floor(W / 2) * math.floor(H / 2)
        
        # Each 2-cm ball occupies 8 (2x2x2) small slots
        small_slots_occupied_by_n2 = n2 * 8
        
        # Remaining slots are filled with 1-cm radius balls
        n1 = total_small_slots - small_slots_occupied_by_n2
        
        energy = n1 * 1 + n2 * 20
        return energy, n1, n2

    # Define the search space for dimensions (in cm, with 0.5 cm steps)
    # A side length of 30cm is a safe upper bound.
    search_dims = [i * 0.5 for i in range(1, 61)]

    # --- Part 1: Box Container Search ---
    for l in search_dims:
        for w in search_dims:
            if w < l: continue  # Ensures l <= w to avoid duplicate shapes
            for h in search_dims:
                if h < w: continue  # Ensures w <= h

                surface_area = 2 * (l * w + l * h + w * h)
                
                if surface_area > max_surface_area:
                    # Since h is the innermost loop and increasing, we can break
                    # if the smallest dimension l is part of a cube-like check.
                    if l == w and w == h:
                        break
                    continue

                energy, n1, n2 = calculate_energy(l, w, h)
                
                if energy > max_energy:
                    max_energy = energy
                    best_config = {
                        "type": "box", "l": l, "w": w, "h": h,
                        "n1": n1, "n2": n2
                    }
            if l == w and best_config and best_config['type'] == 'box' and best_config['l'] == l and best_config['w'] == w and best_config['h'] > max_surface_area / (4*l+2*l):
                break # Optimization for square base prisms

    # --- Part 2: Cylinder Container Search ---
    for r in search_dims:
        for h in search_dims:
            surface_area = 2 * math.pi * r * (r + h)
            if surface_area > max_surface_area:
                continue

            # For a cylinder, we pack within the largest inscribed box.
            # The inscribed box has a square base of side s = r * sqrt(2).
            s = r * math.sqrt(2)
            
            energy, n1, n2 = calculate_energy(s, s, h)

            if energy > max_energy:
                max_energy = energy
                best_config = {
                    "type": "cylinder", "r": r, "h": h,
                    "n1": n1, "n2": n2
                }

    # --- Part 3: Format and Print the Final Answer ---
    if best_config is None:
        print("<<<[0]>>>")
        return

    def format_dim(d):
        """Formats a dimension to int if it's a whole number, else float."""
        return int(d) if d == int(d) else d

    if best_config["type"] == "box":
        # The optimal box found has dimensions l, w, h.
        # We need to check permutations to find the one that was originally tested.
        # However, the calculation is symmetric, so we can use any order.
        # Let's find the permutation that led to the result.
        # The search was l<=w<=h, so we can use that.
        l, w, h = best_config['l'], best_config['w'], best_config['h']
        
        # To be robust, let's re-calculate with all permutations of the best dimensions
        # to find the one that gives the max energy (due to floor function).
        import itertools
        dims = [l, w, h]
        perm_energy = -1
        perm_best = None
        for p in set(itertools.permutations(dims)):
            e, n1_p, n2_p = calculate_energy(p[0], p[1], p[2])
            if e > perm_energy:
                perm_energy = e
                perm_best = {"l": p[0], "w": p[1], "h": p[2], "n1": n1_p, "n2": n2_p}
        
        best_config = perm_best
        
        desc = f"box {format_dim(best_config['l'])}x{format_dim(best_config['w'])}x{format_dim(best_config['h'])}"
        a = best_config["n1"]
        b = best_config["n2"]
        final_answer = f"[{desc}]{a};{b}"
        
    elif best_config["type"] == "cylinder":
        desc = f"cylinder r={format_dim(best_config['r'])}, h={format_dim(best_config['h'])}"
        a = best_config["n1"]
        b = best_config["n2"]
        final_answer = f"[{desc}]{a};{b}"

    print(f"<<<{final_answer}>>>")

solve_pioneer_probe_packing()