import numpy as np
import math

memo = {}

def calculate_max_energy_box(l, w, h):
    """
    Calculates the maximum energy packable in a box of given dimensions
    using a greedy packing strategy. Dimensions are in cm.
    A 2cm-radius ball needs a 4x4x4 cm cube.
    A 1cm-radius ball needs a 2x2x2 cm cube.
    """
    # Use a tuple of sorted dimensions as the key for memoization
    dims = tuple(sorted((l, w, h)))
    if dims in memo:
        return memo[dims]

    # Convert dimensions to 0.5 cm units for the grid
    # A 4cm cube becomes 8 units, a 2cm cube becomes 4 units.
    L, W, H = int(l * 2), int(w * 2), int(h * 2)

    # Check if dimensions are too small even for the smallest ball
    if min(L, W, H) < 4:
        return 0, 0, 0

    grid = np.zeros((L, W, H), dtype=bool)
    
    n2 = 0
    # Pack 2-cm radius balls (8x8x8 grid units) first due to higher energy density
    if min(L, W, H) >= 8:
        for i in range(L - 7):
            for j in range(W - 7):
                for k in range(H - 7):
                    if not np.any(grid[i:i+8, j:j+8, k:k+8]):
                        grid[i:i+8, j:j+8, k:k+8] = True
                        n2 += 1

    n1 = 0
    # Pack 1-cm radius balls (4x4x4 grid units) in the remaining space
    for i in range(L - 3):
        for j in range(W - 3):
            for k in range(H - 3):
                if not np.any(grid[i:i+4, j:j+4, k:k+4]):
                    grid[i:i+4, j:j+4, k:k+4] = True
                    n1 += 1
    
    energy = 20 * n2 + 1 * n1
    result = (energy, n1, n2)
    memo[dims] = result
    return result

def find_optimal_container():
    """
    Searches for the optimal container by iterating through dimensions.
    """
    max_energy = 0
    best_config = None
    
    # Search box containers, WLOG l <= w <= h
    # Dimensions are in 0.5cm steps. Search up to 30cm (60 steps).
    limit = 60
    for i in range(1, limit + 1):
        l = i * 0.5
        for j in range(i, limit + 1):
            w = j * 0.5
            if 2 * (l * w) > 1050: break
            for k in range(j, limit + 1):
                h = k * 0.5
                surface_area = 2 * (l*w + w*h + h*l)
                if surface_area > 1050: break
                
                energy, n1, n2 = calculate_max_energy_box(l, w, h)
                
                if energy > max_energy:
                    max_energy = energy
                    best_config = {"type": "box", "l": l, "w": w, "h": h, "n1": n1, "n2": n2}

    # Search cylinder containers
    # Radius search limit up to 15cm (30 steps)
    r_limit = 30
    for i in range(1, r_limit + 1):
        r = i * 0.5
        if 2 * math.pi * r**2 > 1050: continue
        h_max = (1050 - 2 * math.pi * r**2) / (2 * math.pi * r)
        h_limit = int(h_max / 0.5)

        for j in range(1, h_limit + 1):
            h = j * 0.5
            # Pack into the largest inscribed box as a simplification
            box_side_float = r * math.sqrt(2)
            box_side = math.floor(box_side_float / 0.5) * 0.5
            
            if box_side < 2.0: continue

            energy, n1, n2 = calculate_max_energy_box(box_side, box_side, h)
            
            if energy > max_energy:
                max_energy = energy
                best_config = {"type": "cylinder", "r": r, "h": h, "n1": n1, "n2": n2}

    # Format and print the final answer
    if best_config:
        if best_config["type"] == "box":
            # Ensure consistent order for description
            dims = sorted([best_config['l'], best_config['w'], best_config['h']])
            desc = f"box {dims[0]}x{dims[1]}x{dims[2]}"
        else:
            desc = f"cylinder r={best_config['r']}, h={best_config['h']}"
        
        a = best_config['n1'] # Number of 1-cm balls
        b = best_config['n2'] # Number of 2-cm balls
        
        print(f"[{desc}]{a};{b}")
    else:
        print("[0]")

if __name__ == '__main__':
    find_optimal_container()