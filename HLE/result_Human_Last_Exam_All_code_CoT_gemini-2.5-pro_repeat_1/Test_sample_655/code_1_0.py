import math
import numpy as np

def solve_pioneer_probe_packing():
    """
    This function solves the probe packing problem by:
    1. Defining the three best container shapes based on maximum volume for a given surface area.
    2. Implementing a greedy packing algorithm to fill each container with energy balls.
    3. Comparing the total energy for each container and printing the best result.
    """

    # Step 1: Define the optimal container dimensions found by maximizing volume for A <= 1050 cm^2.
    # Box: 14.0x13.0x12.5 cm (A=1039, V=2275)
    # Cylinder: r=7.0, h=16.5 cm (A=1033.6, V=2539.8)
    # Sphere: r=9.0 cm (A=1017.9, V=3053.6)
    containers = {
        "box": {'dims': (14.0, 13.0, 12.5), 'type': 'box'},
        "cylinder": {'dims': (7.0, 16.5), 'type': 'cylinder'}, # r, h
        "sphere": {'dims': (9.0,), 'type': 'sphere'} # r
    }

    # Define energy balls (large balls are prioritized by placing them first in the list)
    balls_to_pack = [
        {'radius': 2.0, 'energy': 10},
        {'radius': 1.0, 'energy': 1}
    ]
    
    step = 0.5
    max_energy = -1
    best_config = {}

    def is_valid_center(p, r_ball, container):
        """Check if a ball is fully inside the container."""
        x, y, z = p
        ctype = container['type']
        dims = container['dims']
        
        if ctype == 'box':
            L, W, H = dims
            # Box aligned with axes, from 0 to L, 0 to W, 0 to H
            return (r_ball <= x <= L - r_ball) and \
                   (r_ball <= y <= W - r_ball) and \
                   (r_ball <= z <= H - r_ball)
        elif ctype == 'cylinder':
            R_cyl, H_cyl = dims
            # Cylinder axis along z, base centered at (R_cyl, R_cyl)
            return (r_ball <= z <= H_cyl - r_ball) and \
                   ((x - R_cyl)**2 + (y - R_cyl)**2 <= (R_cyl - r_ball)**2)
        elif ctype == 'sphere':
            R_sph = dims[0]
            # Sphere centered at (R_sph, R_sph, R_sph)
            return ((x - R_sph)**2 + (y - R_sph)**2 + (z - R_sph)**2 <= (R_sph - r_ball)**2)
        return False

    def get_search_bounds(container):
        """Get the bounding box for iterating potential centers."""
        ctype = container['type']
        dims = container['dims']
        if ctype == 'box':
            L, W, H = dims
            return (0, L), (0, W), (0, H)
        elif ctype == 'cylinder':
            R, H = dims
            return (0, 2*R), (0, 2*R), (0, H)
        elif ctype == 'sphere':
            R = dims[0]
            return (0, 2*R), (0, 2*R), (0, 2*R)

    # Step 2 & 3: Run greedy packing simulation for each container
    for name, container in containers.items():
        placed_balls = []
        (x_range, y_range, z_range) = get_search_bounds(container)
        
        # Greedy packing algorithm
        for ball_type in balls_to_pack:
            r_ball = ball_type['radius']
            
            ix_coords = np.arange(x_range[0], x_range[1] + step, step)
            iy_coords = np.arange(y_range[0], y_range[1] + step, step)
            iz_coords = np.arange(z_range[0], z_range[1] + step, step)

            for x in ix_coords:
                for y in iy_coords:
                    for z in iz_coords:
                        p = (x, y, z)
                        
                        if not is_valid_center(p, r_ball, container):
                            continue
                            
                        is_overlapping = False
                        for placed_p, placed_r in placed_balls:
                            dist_sq = (p[0] - placed_p[0])**2 + (p[1] - placed_p[1])**2 + (p[2] - placed_p[2])**2
                            if dist_sq < (r_ball + placed_r)**2:
                                is_overlapping = True
                                break
                        
                        if not is_overlapping:
                            placed_balls.append((p, r_ball))

        # Calculate total energy for this container
        n_large = sum(1 for _, r in placed_balls if r == 2.0)
        n_small = sum(1 for _, r in placed_balls if r == 1.0)
        total_energy = 10 * n_large + 1 * n_small

        if total_energy > max_energy:
            max_energy = total_energy
            
            if name == "box":
                desc = f"box {container['dims'][0]}x{container['dims'][1]}x{container['dims'][2]}"
            elif name == "cylinder":
                desc = f"cylinder r={container['dims'][0]}, h={container['dims'][1]}"
            elif name == "sphere":
                desc = f"sphere r={container['dims'][0]}"

            best_config = {
                'desc': desc,
                'n_small': n_small,
                'n_large': n_large
            }

    # Print the final result in the required format
    final_answer_str = f"[{best_config['desc']}]{best_config['n_small']};{best_config['n_large']}"
    print(final_answer_str)

solve_pioneer_probe_packing()