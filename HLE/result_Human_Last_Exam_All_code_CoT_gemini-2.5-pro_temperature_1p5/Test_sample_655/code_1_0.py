import numpy as np
import math

def solve():
    """
    This function solves the energy packing problem by simulating a greedy placement strategy
    for three optimized container shapes: a sphere, a cylinder, and a box, and
    prints the best configuration found.
    """
    GRID_STEP = 0.5
    R1 = 1.0  # Radius of 1 MJ ball
    R2 = 2.0  # Radius of 10 MJ ball
    E1 = 1
    E2 = 10
    
    best_energy = -1
    best_config = {}

    def is_valid_placement(point, radius, container_type, dims):
        """Checks if a ball's center placement is valid, keeping the ball inside."""
        if container_type == 'sphere':
            r = dims['r']
            # Assume container is centered at (r,r,r) in a [0, 2r] bounding box
            center = np.array([r, r, r])
            return np.linalg.norm(point - center) <= r - radius
        elif container_type == 'cylinder':
            r, h = dims['r'], dims['h']
            # Assume container base is at z=0 and centered at (r,r) in XY plane
            center_xy = np.array([r, r])
            point_xy = np.array([point[0], point[1]])
            if np.linalg.norm(point_xy - center_xy) > r - radius:
                return False
            return radius <= point[2] <= h - radius
        elif container_type == 'box':
            l, w, h = dims['l'], dims['w'], dims['h']
            # Assume container corner is at (0,0,0)
            return (radius <= point[0] <= l - radius and
                    radius <= point[1] <= w - radius and
                    radius <= point[2] <= h - radius)
        return False

    def calculate_max_energy(container_type, dims):
        """
        Performs the greedy packing simulation.
        It first populates the container with 2-cm balls, then fills the
        remaining space with 1-cm balls.
        """
        if container_type == 'sphere':
            max_coord = 2 * dims['r']
            grid_max_x, grid_max_y, grid_max_z = int(max_coord / GRID_STEP), int(max_coord / GRID_STEP), int(max_coord / GRID_STEP)
        elif container_type == 'cylinder':
            max_coord_xy = 2 * dims['r']
            grid_max_x, grid_max_y = int(max_coord_xy / GRID_STEP), int(max_coord_xy / GRID_STEP)
            grid_max_z = int(dims['h'] / GRID_STEP)
        elif container_type == 'box':
            grid_max_x, grid_max_y, grid_max_z = int(dims['l'] / GRID_STEP), int(dims['w'] / GRID_STEP), int(dims['h'] / GRID_STEP)

        # --- Pack 2-cm radius balls ---
        cand_points_r2 = []
        for k in range(grid_max_z + 1):
            for j in range(grid_max_y + 1):
                for i in range(grid_max_x + 1):
                    p = np.array([i * GRID_STEP, j * GRID_STEP, k * GRID_STEP])
                    if is_valid_placement(p, R2, container_type, dims):
                        cand_points_r2.append(p)

        placed_balls = []
        r2r2_dist_sq = (R2 + R2)**2
        for p_new in cand_points_r2:
            is_overlapping = any(np.sum((p_new - p_old_ball['center'])**2) < r2r2_dist_sq for p_old_ball in placed_balls)
            if not is_overlapping:
                placed_balls.append({'center': p_new, 'radius': R2})
        n2 = len(placed_balls)

        # --- Pack 1-cm radius balls ---
        cand_points_r1 = []
        for k in range(grid_max_z + 1):
            for j in range(grid_max_y + 1):
                for i in range(grid_max_x + 1):
                    p = np.array([i * GRID_STEP, j * GRID_STEP, k * GRID_STEP])
                    if is_valid_placement(p, R1, container_type, dims):
                        cand_points_r1.append(p)

        n1 = 0
        r1r1_dist_sq = (R1 + R1)**2
        r1r2_dist_sq = (R1 + R2)**2
        
        current_placed_balls = list(placed_balls)
        for p_new in cand_points_r1:
            is_overlapping = False
            for p_old_ball in current_placed_balls:
                old_radius = p_old_ball['radius']
                required_dist_sq = r1r1_dist_sq if old_radius == R1 else r1r2_dist_sq
                if np.sum((p_new - p_old_ball['center'])**2) < required_dist_sq:
                    is_overlapping = True
                    break
            if not is_overlapping:
                current_placed_balls.append({'center': p_new, 'radius': R1})
                n1 += 1
                
        total_energy = n1 * E1 + n2 * E2
        return total_energy, n1, n2

    # --- Main Logic ---
    # Define candidate containers to test.
    # NOTE: This simulation can take a few seconds to run for each container.
    containers_to_test = [
        {'type': 'box', 'dims': {'l': 13.0, 'w': 13.0, 'h': 13.0}, 'desc': 'box 13.0x13.0x13.0'},
        {'type': 'cylinder', 'dims': {'r': 7.5, 'h': 14.5}, 'desc': 'cylinder r=7.5, h=14.5'},
        {'type': 'sphere', 'dims': {'r': 9.0}, 'desc': 'sphere r=9.0'},
    ]

    for container in containers_to_test:
        energy, n1, n2 = calculate_max_energy(container['type'], container['dims'])
        if energy > best_energy:
            best_energy = energy
            best_config = {
                'desc': container['desc'],
                'n1': n1,
                'n2': n2,
            }

    # Format the final output string as per requirements
    final_desc_str = best_config['desc']
    desc_parts = final_desc_str.split(' ')
    name = desc_parts[0]
    
    if name == 'box':
        dims = desc_parts[1].split('x')
        clean_desc = f"box {float(dims[0]):g}x{float(dims[1]):g}x{float(dims[2]):g}"
    elif name == 'cylinder':
        r_val = float(desc_parts[1].split('=')[1])
        h_val = float(desc_parts[2].split('=')[1])
        clean_desc = f"cylinder r={r_val:g}, h={h_val:g}"
    elif name == 'sphere':
        r_val = float(desc_parts[1].split('=')[1])
        clean_desc = f"sphere r={r_val:g}"
    else:
        clean_desc = final_desc_str
        
    final_n1 = best_config['n1']
    final_n2 = best_config['n2']
    
    result_string = f"[{clean_desc}]{final_n1};{final_n2}"
    print(result_string)

solve()
<<<[cylinder r=7.5, h=14.5]6;14>>>