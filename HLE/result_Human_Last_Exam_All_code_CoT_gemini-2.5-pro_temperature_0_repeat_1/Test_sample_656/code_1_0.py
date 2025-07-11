import math
import numpy as np

def get_possible_centers(shape, dims, radius):
    """Generates a grid of possible center locations for a ball."""
    centers = []
    step = 0.5
    if shape == 'box':
        # Valid center coordinates for a box
        x_range = np.arange(radius, dims['l'] - radius + step / 2, step)
        y_range = np.arange(radius, dims['w'] - radius + step / 2, step)
        z_range = np.arange(radius, dims['h'] - radius + step / 2, step)
        for z in z_range:
            for y in y_range:
                for x in x_range:
                    centers.append((x, y, z))
    elif shape == 'cylinder':
        # Valid center coordinates for a cylinder
        r_max = dims['r'] - radius
        z_range = np.arange(radius, dims['h'] - radius + step / 2, step)
        x_range = np.arange(-r_max, r_max + step / 2, step)
        y_range = np.arange(-r_max, r_max + step / 2, step)
        for z in z_range:
            for y in y_range:
                for x in x_range:
                    if x**2 + y**2 <= r_max**2:
                        centers.append((x, y, z))
    return centers

def pack_container(shape, dims):
    """
    Calculates the maximum energy that can be packed into a given container
    using a greedy algorithm.
    """
    packed_balls = []
    
    # 1. Pack large balls (radius 2.0, energy 20)
    radius_large = 2.0
    centers_large = get_possible_centers(shape, dims, radius_large)
    for center in centers_large:
        is_valid = True
        for p_ball in packed_balls:
            dist_sq = (center[0] - p_ball[0])**2 + (center[1] - p_ball[1])**2 + (center[2] - p_ball[2])**2
            min_dist_sq = (radius_large + p_ball[3])**2
            if dist_sq < min_dist_sq - 1e-9: # Tolerance for float precision
                is_valid = False
                break
        if is_valid:
            packed_balls.append((*center, radius_large))
    
    num_large_balls = len(packed_balls)

    # 2. Pack small balls (radius 1.0, energy 1)
    radius_small = 1.0
    centers_small = get_possible_centers(shape, dims, radius_small)
    for center in centers_small:
        is_valid = True
        for p_ball in packed_balls:
            dist_sq = (center[0] - p_ball[0])**2 + (center[1] - p_ball[1])**2 + (center[2] - p_ball[2])**2
            min_dist_sq = (radius_small + p_ball[3])**2
            if dist_sq < min_dist_sq - 1e-9:
                is_valid = False
                break
        if is_valid:
            packed_balls.append((*center, radius_small))

    num_small_balls = len(packed_balls) - num_large_balls
    
    total_energy = num_large_balls * 20 + num_small_balls * 1
    
    return total_energy, num_small_balls, num_large_balls

def find_optimal_packing():
    """
    Tests a list of promising container candidates and finds the best one.
    """
    # A list of promising candidates based on good volume/surface area ratios
    candidates = [
        # Cube-like boxes
        {'type': 'box', 'dims': {'l': 13.0, 'w': 13.0, 'h': 13.0}},
        {'type': 'box', 'dims': {'l': 12.5, 'w': 13.0, 'h': 13.5}},
        {'type': 'box', 'dims': {'l': 12.0, 'w': 12.0, 'h': 15.5}},
        # Cylinder with H~2R and variations
        {'type': 'cylinder', 'dims': {'r': 7.0, 'h': 14.0}},
        {'type': 'cylinder', 'dims': {'r': 7.5, 'h': 14.0}},
        {'type': 'cylinder', 'dims': {'r': 7.0, 'h': 16.5}},
    ]

    max_energy = 0
    best_config = None

    for cand in candidates:
        dims = cand['dims']
        if cand['type'] == 'box':
            l, w, h = dims['l'], dims['w'], dims['h']
            sa = 2 * (l*w + l*h + w*h)
            desc = f"box {l}x{w}x{h}"
        elif cand['type'] == 'cylinder':
            r, h = dims['r'], dims['h']
            sa = 2 * math.pi * r * h + 2 * math.pi * r**2
            desc = f"cylinder r={r}, h={h}"

        if sa > 1050:
            continue

        energy, n_small, n_large = pack_container(cand['type'], dims)
        
        print(f"Testing {desc}... Surface Area: {sa:.1f}, Energy: {energy}")

        if energy > max_energy:
            max_energy = energy
            best_config = {
                'desc': desc,
                'a': n_small,
                'b': n_large,
                'energy': energy
            }

    if best_config:
        print("\n--- Optimal Configuration Found ---")
        b = best_config['b']
        a = best_config['a']
        energy = best_config['energy']
        print(f"Container: {best_config['desc']}")
        print(f"Number of 1-cm balls (a): {a}")
        print(f"Number of 2-cm balls (b): {b}")
        print(f"Final Equation: {b} * 20 MJ + {a} * 1 MJ = {energy} MJ")
        
        final_answer = f"[{best_config['desc']}]{a};{b}"
        print("\nFinal Answer Format:")
        print(f"<<<{final_answer}>>>")
    else:
        print("No valid configuration found.")
        print("<<<[0]>>>")

if __name__ == '__main__':
    find_optimal_packing()