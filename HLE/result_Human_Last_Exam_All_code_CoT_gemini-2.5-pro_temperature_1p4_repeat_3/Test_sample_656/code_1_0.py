import numpy as np
import math

def solve_packing_problem():
    """
    This function solves the container design problem by performing a computational search.
    It explores different container shapes and dimensions to find the one that can hold
    the maximum amount of energy, respecting the given constraints.
    """

    # --- Geometric and Packing Helper Functions ---
    def check_overlap(center1, radius1, center2, radius2):
        """Checks if two spheres overlap by comparing the squared distance between their centers to the squared sum of their radii."""
        dist_sq = (center1[0] - center2[0])**2 + (center1[1] - center2[1])**2 + (center1[2] - center2[2])**2
        return dist_sq < (radius1 + radius2)**2

    def is_inside_box(center, radius, length, width, height):
        """Checks if a sphere is fully contained within a box defined in the first octant."""
        cx, cy, cz = center
        return (radius <= cx <= length - radius and
                radius <= cy <= width - radius and
                radius <= cz <= height - radius)

    def is_inside_cylinder(center, radius, cyl_radius, cyl_height):
        """Checks if a sphere is fully contained within a cylinder. The cylinder's base is centered at (R,R) in the xy-plane."""
        cx, cy, cz = center
        if not (radius <= cz <= cyl_height - radius):
            return False
        if math.sqrt((cx - cyl_radius)**2 + (cy - cyl_radius)**2) + radius > cyl_radius:
            return False
        return True

    def pack_container(container_type, dims):
        """
        Simulates packing energy balls into a given container using a greedy, grid-based approach.
        It prioritizes placing the larger, more valuable 2-cm balls first.
        """
        placed_balls = []

        # Define the grid for ball centers based on container dimensions
        if container_type == 'box':
            length, width, height = dims
            ranges = [np.arange(0.0, d + 0.5, 0.5) for d in dims]
            is_inside_func = lambda c, r: is_inside_box(c, r, length, width, height)
        else:  # 'cylinder'
            cyl_radius, cyl_height = dims
            ranges = [np.arange(0.0, 2 * cyl_radius + 0.5, 0.5)] * 2 + [np.arange(0.0, cyl_height + 0.5, 0.5)]
            is_inside_func = lambda c, r: is_inside_cylinder(c, r, cyl_radius, cyl_height)

        # --- Pass 1: Pack 2-cm radius balls (value 20 MJ) ---
        radius_2cm = 2.0
        for cz in ranges[2]:
            for cy in ranges[1]:
                for cx in ranges[0]:
                    center = (cx, cy, cz)
                    if is_inside_func(center, radius_2cm):
                        # Check for overlap with already placed balls
                        if not any(check_overlap(center, radius_2cm, p_center, p_radius) for p_center, p_radius in placed_balls):
                            placed_balls.append((center, radius_2cm))
        num_2cm_balls = len(placed_balls)

        # --- Pass 2: Pack 1-cm radius balls (value 1 MJ) ---
        radius_1cm = 1.0
        for cz in ranges[2]:
            for cy in ranges[1]:
                for cx in ranges[0]:
                    center = (cx, cy, cz)
                    if is_inside_func(center, radius_1cm):
                        if not any(check_overlap(center, radius_1cm, p_center, p_radius) for p_center, p_radius in placed_balls):
                            placed_balls.append((center, radius_1cm))
                            
        num_1cm_balls = len(placed_balls) - num_2cm_balls
        return num_2cm_balls, num_1cm_balls

    # --- Main Search Loop ---
    max_energy = 0
    best_config = None
    max_sa = 1050.0

    # Search a few promising box dimensions based on preliminary analysis
    box_dims_to_check = [
        (16.0, 16.0, 8.0), 
        (28.0, 8.0, 8.0)
    ]
    for l, w, h in set(box_dims_to_check):
        sa = 2 * (l*w + l*h + w*h)
        if sa > max_sa:
            continue
        num_b, num_a = pack_container('box', (l, w, h))
        energy = 20 * num_b + num_a
        if energy > max_energy:
            max_energy = energy
            best_config = {"type": "box", "l": l, "w": w, "h": h, "a": num_a, "b": num_b}

    # Search promising cylinder dimensions. Analysis suggests an optimal shape is when height is twice the radius.
    cyl_dims_to_check = [
        (7.0, 14.0),
        (7.5, 13.0) # SA for r=7.5,h=15 is >1050
    ]
    for r, h in cyl_dims_to_check:
        sa = 2 * math.pi * r * (r + h)
        if sa > max_sa:
            continue
        num_b, num_a = pack_container('cylinder', (r, h))
        energy = 20 * num_b + num_a
        if energy > max_energy:
            max_energy = energy
            best_config = {"type": "cylinder", "r": r, "h": h, "a": num_a, "b": num_b}
    
    # --- Format and Print Final Answer ---
    if best_config:
        cfg = best_config
        if cfg["type"] == "box":
            desc = f'box {cfg["l"]:.1f}x{cfg["w"]:.1f}x{cfg["h"]:.1f}'
            print(f"[{desc}]{cfg['a']};{cfg['b']}")
        else: # cylinder
            desc = f'cylinder r={cfg["r"]:.1f}, h={cfg["h"]:.1f}'
            print(f"[{desc}]{cfg['a']};{cfg['b']}")
    else:
        print("[0]")

solve_packing_problem()