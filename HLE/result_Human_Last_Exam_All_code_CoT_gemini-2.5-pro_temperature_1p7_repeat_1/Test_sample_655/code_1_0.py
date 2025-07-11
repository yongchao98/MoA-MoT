import math

def solve_packing_problem():
    """
    This function designs the optimal container and packs energy balls to maximize energy.
    It focuses on the most promising container (a sphere) and uses a greedy packing algorithm.
    """
    
    # Based on preliminary analysis, the sphere is the most volume-efficient container.
    # SA = 4*pi*r^2 <= 1050  => r <= 9.14. Max r that is multiple of 0.5 is r=9.0
    # SA for r=9.0 is 4*pi*81 ~= 1017.88, which is valid.
    container = {
        'type': 'sphere', 
        'dims': [9.0], 
        'desc': 'sphere r=9.0'
    }

    # Sphere is centered at (0,0,0)
    R = container['dims'][0]
    min_coord = -R
    max_coord = R
        
    placed_balls = []

    # --- Step 1: Pack 2-cm radius balls (high-energy) ---
    r_large = 2.0
    
    # Generate a list of all valid candidate center points for large balls
    candidates_large = []
    # Iterate through the grid points inside the container's bounding box
    for i in range(int(2 * min_coord), int(2 * max_coord) + 1):
        for j in range(int(2 * min_coord), int(2 * max_coord) + 1):
            for k in range(int(2 * min_coord), int(2 * max_coord) + 1):
                center = (i * 0.5, j * 0.5, k * 0.5)
                # Check if this center is valid for a large ball (sphere centered at origin)
                if math.sqrt(center[0]**2 + center[1]**2 + center[2]**2) + r_large <= R:
                    candidates_large.append(center)
    
    # To get a denser packing, sort candidates from the center outwards
    candidates_large.sort(key=lambda p: p[0]**2 + p[1]**2 + p[2]**2)

    # Greedy placement for large balls
    min_dist_sq_large = (r_large + r_large)**2
    for p_new in candidates_large:
        can_place = True
        for p_placed in placed_balls:
            dist_sq = (p_new[0] - p_placed[0])**2 + (p_new[1] - p_placed[1])**2 + (p_new[2] - p_placed[2])**2
            if dist_sq < min_dist_sq_large:
                can_place = False
                break
        if can_place:
            placed_balls.append(p_new)

    num_2cm_balls = len(placed_balls)
    # Convert placed centers to a list of (x,y,z,r) tuples for the next phase
    placed_balls_with_radius = [ (p[0], p[1], p[2], r_large) for p in placed_balls ]

    # --- Step 2: Pack 1-cm radius balls (low-energy) in the gaps ---
    r_small = 1.0

    # Generate candidate centers for small balls
    candidates_small = []
    for i in range(int(2 * min_coord), int(2 * max_coord) + 1):
        for j in range(int(2 * min_coord), int(2 * max_coord) + 1):
            for k in range(int(2 * min_coord), int(2 * max_coord) + 1):
                center = (i * 0.5, j * 0.5, k * 0.5)
                if math.sqrt(center[0]**2 + center[1]**2 + center[2]**2) + r_small <= R:
                    candidates_small.append(center)

    # Sort candidates from center outwards
    candidates_small.sort(key=lambda p: p[0]**2 + p[1]**2 + p[2]**2)

    # Greedy placement for small balls, checking against all previously placed balls
    for p_new in candidates_small:
        can_place = True
        for p_placed in placed_balls_with_radius:
            r_placed = p_placed[3]
            min_dist_sq = (r_small + r_placed)**2
            dist_sq = (p_new[0] - p_placed[0])**2 + (p_new[1] - p_placed[1])**2 + (p_new[2] - p_placed[2])**2
            if dist_sq < min_dist_sq:
                can_place = False
                break
        if can_place:
            placed_balls_with_radius.append((p_new[0], p_new[1], p_new[2], r_small))

    num_1cm_balls = len(placed_balls_with_radius) - num_2cm_balls
    
    # Final output
    container_desc = container['desc']
    a = num_1cm_balls
    b = num_2cm_balls
    
    print(f"[{container_desc}]{a};{b}")

if __name__ == '__main__':
    solve_packing_problem()
<<<[sphere r=9.0]144;57>>>