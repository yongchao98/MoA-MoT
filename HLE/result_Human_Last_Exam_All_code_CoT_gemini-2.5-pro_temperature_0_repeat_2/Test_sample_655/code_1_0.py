import math

def solve_packing_problem():
    """
    This function calculates the optimal packing of energy balls in a spherical
    container using a greedy algorithm.
    """
    # 1. Define the optimal container: a sphere with r=9.0 cm.
    container_radius = 9.0
    container_desc = f"sphere r={container_radius}"
    grid_step = 0.5

    # 2. Generate potential center points for the balls.
    # The points are on a 0.5cm grid and sorted by distance from the origin
    # to prioritize packing from the center outwards.
    potential_centers = []
    max_coord_idx = int(container_radius / grid_step)

    for i in range(-max_coord_idx, max_coord_idx + 1):
        for j in range(-max_coord_idx, max_coord_idx + 1):
            for k in range(-max_coord_idx, max_coord_idx + 1):
                pos = (i * grid_step, j * grid_step, k * grid_step)
                dist_sq = pos[0]**2 + pos[1]**2 + pos[2]**2
                # Add point if it's within the container's bounding box
                if dist_sq <= container_radius**2:
                    potential_centers.append({'pos': pos, 'dist_sq': dist_sq})

    potential_centers.sort(key=lambda p: p['dist_sq'])

    # Helper function to check if a ball can be placed at a position
    def is_valid_placement(center, radius, placed_balls):
        # Check if the ball is fully inside the container
        if math.sqrt(center[0]**2 + center[1]**2 + center[2]**2) + radius > container_radius:
            return False
        # Check for overlap with already placed balls
        for placed in placed_balls:
            dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(center, placed['center']))
            min_dist = radius + placed['radius']
            # Use a small tolerance for floating point comparisons
            if dist_sq < min_dist**2 - 1e-9:
                return False
        return True

    # 3. Run the greedy packing algorithm.
    placed_balls = []
    
    # Phase 1: Pack 2-cm radius balls (10 MJ)
    radius_2cm = 2.0
    for p in potential_centers:
        center = p['pos']
        if is_valid_placement(center, radius_2cm, placed_balls):
            placed_balls.append({'center': center, 'radius': radius_2cm})
    
    n2 = len(placed_balls)

    # Phase 2: Pack 1-cm radius balls (1 MJ) in the remaining space
    radius_1cm = 1.0
    # We can check all potential centers again, the overlap check will handle it
    for p in potential_centers:
        center = p['pos']
        if is_valid_placement(center, radius_1cm, placed_balls):
            placed_balls.append({'center': center, 'radius': radius_1cm})

    n1 = len(placed_balls) - n2
    
    # 4. Output the results
    total_energy = 10 * n2 + 1 * n1
    print(f"Calculation for container: {container_desc}")
    print(f"Number of 2-cm balls (b): {n2}")
    print(f"Number of 1-cm balls (a): {n1}")
    print(f"The final equation for energy is: 10 * {n2} + 1 * {n1} = {total_energy} MJ")
    print("\nFinal Answer String:")
    print(f"[{container_desc}]{n1};{n2}")


solve_packing_problem()