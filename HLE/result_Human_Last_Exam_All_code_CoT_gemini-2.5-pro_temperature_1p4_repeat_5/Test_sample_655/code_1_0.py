import math
import numpy as np

# Global tracker for the best found configuration
best_config = {
    "energy": -1,
    "description": "",
    "n_small": 0,
    "n_big": 0
}

def distance_sq(p1, p2):
    """Calculates the squared Euclidean distance between two 3D points."""
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

def pack_container(shape, dims):
    """
    Performs a greedy packing simulation for a given container.
    It packs larger, more energetic balls first.
    It prioritizes placing balls towards the geometric center of the container.
    """
    placed_balls = []
    
    # 1. Define container-specific properties and boundary checks
    is_in_bounds_func = None
    container_center = None
    bounding_box_ranges = ()

    if shape == 'box':
        l, w, h = dims['l'], dims['w'], dims['h']
        def is_in_bounds(center, radius):
            cx, cy, cz = center
            return (radius <= cx <= l - radius and
                    radius <= cy <= w - radius and
                    radius <= cz <= h - radius)
        is_in_bounds_func = is_in_bounds
        container_center = (l / 2.0, w / 2.0, h / 2.0)
        bounding_box_ranges = (np.arange(0, l + 0.5, 0.5),
                               np.arange(0, w + 0.5, 0.5),
                               np.arange(0, h + 0.5, 0.5))
    elif shape == 'cylinder':
        r, h = dims['r'], dims['h']
        # The cylinder is placed with its base on the z=0 plane,
        # centered at (r, r) in the xy plane.
        def is_in_bounds(center, radius):
            cx, cy, cz = center
            # Using a small tolerance for floating point comparisons
            return (math.sqrt((cx - r)**2 + (cy - r)**2) + radius <= r + 1e-9 and
                    radius <= cz <= h - radius)
        is_in_bounds_func = is_in_bounds
        container_center = (r, r, h / 2.0)
        bounding_box_ranges = (np.arange(0, 2 * r + 0.5, 0.5),
                               np.arange(0, 2 * r + 0.5, 0.5),
                               np.arange(0, h + 0.5, 0.5))
    elif shape == 'sphere':
        r = dims['r']
        # The sphere is centered at (r, r, r).
        def is_in_bounds(center, radius):
            cx, cy, cz = center
            return math.sqrt((cx - r)**2 + (cy - r)**2 + (cz - r)**2) + radius <= r + 1e-9
        is_in_bounds_func = is_in_bounds
        container_center = (r, r, r)
        bounding_box_ranges = (np.arange(0, 2 * r + 0.5, 0.5),
                               np.arange(0, 2 * r + 0.5, 0.5),
                               np.arange(0, 2 * r + 0.5, 0.5))

    # 2. Create and sort a list of potential center points for the balls
    potential_centers = []
    for cx in bounding_box_ranges[0]:
        for cy in bounding_box_ranges[1]:
            for cz in bounding_box_ranges[2]:
                # Pre-filter points that are obviously outside the shape
                if is_in_bounds_func((cx,cy,cz), 0):
                    potential_centers.append((cx, cy, cz))
    
    potential_centers.sort(key=lambda p: distance_sq(p, container_center))

    # 3. Greedy packing loop (big balls first, then small balls)
    ball_types = [(2, 10), (1, 1)]  # (radius, energy)
    n_big, n_small = 0, 0

    for radius, _ in ball_types:
        temp_placed_count = 0
        for center in potential_centers:
            if not is_in_bounds_func(center, radius):
                continue

            is_overlapping = False
            for ball in placed_balls:
                required_dist_sq = (radius + ball['radius'])**2
                if distance_sq(center, ball['center']) < required_dist_sq - 1e-9:
                    is_overlapping = True
                    break
            
            if not is_overlapping:
                placed_balls.append({'center': center, 'radius': radius})
                temp_placed_count += 1
        
        if radius == 2:
            n_big = temp_placed_count
        else:
            n_small = temp_placed_count

    total_energy = 10 * n_big + 1 * n_small
    # The number of small balls is the total minus the number of big balls
    actual_n_small = len(placed_balls) - n_big
    return total_energy, n_big, actual_n_small

def solve():
    """
    Main function to search for the optimal container.
    """
    global best_config
    step = 0.5
    max_sa = 1050.0
    heuristic_density = (10 / (4/3 * math.pi * 2**3)) * 0.5

    # --- Box Search (l >= w >= h) ---
    h_max = (max_sa / 6)**0.5
    for h in np.arange(4.0, h_max + step, step):
        w_max = -h + math.sqrt(h**2 + max_sa / 2)
        for w in np.arange(h, w_max + step, step):
            if 2 * w * h > max_sa: continue
            l_max = (max_sa / 2 - w * h) / (w + h)
            for l in np.arange(w, l_max + step, step):
                if 2 * (l*w + w*h + h*l) > max_sa: continue
                volume = l * w * h
                if volume * heuristic_density < best_config["energy"]: continue
                dims = {'l': l, 'w': w, 'h': h}
                energy, n_big, n_small = pack_container('box', dims)
                if energy > best_config["energy"]:
                    best_config = {"energy": energy, "description": f"box {l}x{w}x{h}", "n_small": n_small, "n_big": n_big}

    # --- Cylinder Search ---
    r_max = math.sqrt(max_sa / (2 * math.pi))
    for r in np.arange(2.0, r_max + step, step):
        sa_base = 2 * math.pi * r**2
        if sa_base > max_sa: continue
        h_max = (max_sa - sa_base) / (2 * math.pi * r)
        for h in np.arange(4.0, h_max + step, step):
            volume = math.pi * r**2 * h
            if volume * heuristic_density < best_config["energy"]: continue
            dims = {'r': r, 'h': h}
            energy, n_big, n_small = pack_container('cylinder', dims)
            if energy > best_config["energy"]:
                 best_config = {"energy": energy, "description": f"cylinder r={r}, h={h}", "n_small": n_small, "n_big": n_big}

    # --- Sphere Search ---
    r_max = math.sqrt(max_sa / (4 * math.pi))
    for r in np.arange(2.0, r_max + step, step):
        if 4 * math.pi * r**2 > max_sa: continue
        volume = 4/3 * math.pi * r**3
        if volume * heuristic_density < best_config["energy"]: continue
        dims = {'r': r}
        energy, n_big, n_small = pack_container('sphere', dims)
        if energy > best_config["energy"]:
            best_config = {"energy": energy, "description": f"sphere r={r}", "n_small": n_small, "n_big": n_big}
    
    print(f'[{best_config["description"]}]{best_config["n_small"]};{best_config["n_big"]}')

if __name__ == '__main__':
    solve()