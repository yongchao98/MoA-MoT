import math
import numpy as np

def solve_packing_problem():
    """
    Solves the exo-planet probe packing problem by searching for the optimal
    container and ball configuration to maximize energy.
    """

    def floor_to_step(value, step=0.5):
        """Rounds a value down to the nearest multiple of the given step."""
        if value is None or value < step:
            return 0
        return math.floor(value / step) * step

    def pack_balls_heuristic(l, w, h):
        """
        A heuristic to calculate the number of balls that can be packed into a box.
        It prioritizes large balls in a grid and fills remaining space and voids
        with small balls.
        """
        # Ensure dimensions are ordered for consistent slab calculation
        l, w, h = sorted((l, w, h), reverse=True)

        # 1. Pack large balls (radius=2, diameter=4) in a simple grid
        num_l = l // 4
        num_w = w // 4
        num_h = h // 4
        n_large = int(num_l * num_w * num_h)

        # Calculate volume occupied by the large ball grid
        vol_large_l = num_l * 4.0
        vol_large_w = num_w * 4.0
        vol_large_h = num_h * 4.0

        n_small = 0

        # 2. Pack small balls (radius=1, diameter=2) in interstitial voids
        # A small ball (r=1) fits in the central void between 8 large balls (r=2)
        if num_l > 1 and num_w > 1 and num_h > 1:
            n_small_voids = (num_l - 1) * (num_w - 1) * (num_h - 1)
            n_small += n_small_voids

        # 3. Pack small balls in the three leftover slabs of space
        # This simplified model avoids double-counting corners.
        rem_l = l - vol_large_l
        if rem_l >= 2.0:
            # Slab on the side of length l
            n_small += (rem_l // 2) * (w // 2) * (h // 2)

        rem_w = w - vol_large_w
        if rem_w >= 2.0:
            # Slab on the side of width w (but only over the packed large balls area)
            n_small += vol_large_l // 2 * (rem_w // 2) * (h // 2)

        rem_h = h - vol_large_h
        if rem_h >= 2.0:
            # Slab on the side of height h (but only over the packed large balls area)
            n_small += (vol_large_l // 2) * (vol_large_w // 2) * (rem_h // 2)

        energy = n_large * 20 + n_small * 1
        return energy, n_large, int(n_small)

    best_config = {
        'energy': 0, 'shape': None, 'dims': None, 'n_small': 0, 'n_large': 0
    }
    surface_area_limit = 1050
    step = 0.5
    # A reasonable upper bound for any single dimension
    max_dim_check = int(surface_area_limit / (4*step) / 2) # e.g. a flat box 0.5*X*2SA = X

    dims_to_check = np.arange(step, max_dim_check + step, step)
    
    # --- Search Box Shapes (l >= w >= h) ---
    for h in dims_to_check:
        for w in dims_to_check:
            if w < h: continue
            
            numerator = surface_area_limit - 2 * w * h
            denominator = 2 * w + 2 * h
            if numerator <= 0: break

            l_max = numerator / denominator
            if l_max < w: continue
            l = floor_to_step(l_max, step)
            if l == 0: continue
            
            energy, n_large, n_small = pack_balls_heuristic(l, w, h)
            
            if energy > best_config['energy']:
                best_config = {'energy': energy, 'shape': 'box', 'dims': (l, w, h), 
                               'n_small': n_small, 'n_large': n_large}
    
    # --- Search Cylinder Shapes ---
    for r in dims_to_check:
        numerator = surface_area_limit - 2 * math.pi * r**2
        if numerator <= 0: break
            
        denominator = 2 * math.pi * r
        h_max = numerator / denominator
        h = floor_to_step(h_max, step)
        if h == 0: continue
            
        # For cylinders, we approximate packing efficiency by using the largest
        # possible square-base box that fits inside. This is a conservative heuristic.
        box_side = floor_to_step((2*r) / math.sqrt(2), step)
        energy, n_large, n_small = pack_balls_heuristic(box_side, box_side, h)

        if energy > best_config['energy']:
             best_config = {'energy': energy, 'shape': 'cylinder', 'dims': (r, h),
                            'n_small': n_small, 'n_large': n_large}

    # --- Format and Print Final Answer ---
    shape = best_config['shape']
    dims = best_config['dims']
    n_small = best_config['n_small']
    n_large = best_config['n_large']

    description = "No solution found"
    if shape == 'box':
        # Ensure consistent dimension order: longest to shortest
        l, w, h = sorted(dims, reverse=True)
        description = f"{shape} {l}x{w}x{h}"
    elif shape == 'cylinder':
        r, h = dims
        description = f"{shape} r={r}, h={h}"

    final_answer = f"[{description}]{n_small};{n_large}"
    print(final_answer)
    return final_answer

# Run the solver
final_answer_string = solve_packing_problem()
# The required final answer block
print(f"\n<<<{final_answer_string}>>>")