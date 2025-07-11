import math

def get_max_balls(length, width, height):
    """
    Calculates the number of balls that can fit in a box using a greedy algorithm.
    The algorithm iterates through possible center locations on a 0.5cm grid and places a ball
    if it does not collide with any existing balls.
    """
    ball_radius = 2.0
    min_dist_sq = (2 * ball_radius)**2

    # Define the valid range for ball centers. A center must be at least `ball_radius` from any wall.
    # The step for coordinates is 0.5 cm.
    x_coords = [i * 0.5 for i in range(int(ball_radius / 0.5), int((length - ball_radius) / 0.5) + 1)]
    y_coords = [i * 0.5 for i in range(int(ball_radius / 0.5), int((width - ball_radius) / 0.5) + 1)]
    z_coords = [i * 0.5 for i in range(int(ball_radius / 0.5), int(height - ball_radius) / 0.5) + 1)]

    # If any dimension is too small, no balls can be placed.
    if not x_coords or not y_coords or not z_coords:
        return 0

    placed_centers = []
    # Iterate through all possible center points in a systematic order (z, then y, then x).
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                candidate_center = (x, y, z)
                can_place = True
                # Check for collision with already placed balls.
                for placed_center in placed_centers:
                    dist_sq = (candidate_center[0] - placed_center[0])**2 + \
                              (candidate_center[1] - placed_center[1])**2 + \
                              (candidate_center[2] - placed_center[2])**2
                    if dist_sq < min_dist_sq:
                        can_place = False
                        break
                
                if can_place:
                    placed_centers.append(candidate_center)
    
    return len(placed_centers)

def find_optimal_box():
    """
    Searches for an optimal box with integer dimensions that can hold at least 27 balls
    and has a surface area less than the original 12x12x12 box.
    """
    initial_dim = 12
    initial_surface_area = 6 * (initial_dim ** 2)
    min_balls_required = 27

    best_solution = None
    min_surface_area = initial_surface_area

    # Search integer dimensions within a reasonable range.
    # To avoid duplicates, we search with l <= w <= h.
    for l in range(4, 15): 
        for w in range(l, 17): 
            for h in range(w, 19): 
                
                current_surface_area = 2 * (l*w + l*h + w*h)

                if current_surface_area >= min_surface_area:
                    continue
                
                num_balls = get_max_balls(float(l), float(w), float(h))
                
                if num_balls >= min_balls_required:
                    min_surface_area = current_surface_area
                    best_solution = (l, w, h, int(current_surface_area))

    if best_solution:
        l, w, h, area = best_solution
        print(f"{l}:{w}:{h}:{area}")
    else:
        print("0")

# Run the optimization search.
find_optimal_box()