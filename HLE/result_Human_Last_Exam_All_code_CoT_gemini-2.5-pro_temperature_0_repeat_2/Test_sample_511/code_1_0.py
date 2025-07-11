import math

def solve_pioneer_packing():
    """
    Finds the optimal container dimensions to pack energy balls.
    """
    # --- Problem Parameters ---
    ball_radius = 2.0
    min_dist_between_centers = 2 * ball_radius
    min_dist_sq = min_dist_between_centers**2
    precision = 0.5

    # --- Original Box Configuration ---
    orig_dims = (12, 12, 12)
    orig_surface_area = 2 * (orig_dims[0] * orig_dims[1] + 
                             orig_dims[1] * orig_dims[2] + 
                             orig_dims[2] * orig_dims[0])
    
    # Calculate the number of balls in the original box using simple cubic packing
    # This determines the minimum number of balls the new box must hold.
    nx = math.floor((orig_dims[0] - 2 * ball_radius) / min_dist_between_centers) + 1
    ny = math.floor((orig_dims[1] - 2 * ball_radius) / min_dist_between_centers) + 1
    nz = math.floor((orig_dims[2] - 2 * ball_radius) / min_dist_between_centers) + 1
    min_balls_required = nx * ny * nz

    # --- Search for a more efficient box ---
    best_box_dims = None
    min_surface_area = orig_surface_area

    # Define the search space for integer dimensions L, W, H.
    # We search in a reasonable range around the original 12cm dimension.
    # To avoid duplicate permutations (e.g., 10x11x12 and 12x10x11), we enforce L >= W >= H.
    for l in range(8, 16):
        for w in range(8, l + 1):
            for h in range(8, w + 1):
                
                current_surface_area = 2 * (l*w + w*h + h*l)
                
                # Pruning Step 1: If the surface area is not an improvement, skip.
                if current_surface_area >= min_surface_area:
                    continue

                # --- Greedy Packing Simulation ---
                
                # Generate all possible center coordinates based on precision.
                # A ball center must be at least `radius` away from each wall.
                x_coords = [i * precision for i in range(int(ball_radius / precision), int((l - ball_radius) / precision) + 1)]
                y_coords = [i * precision for i in range(int(ball_radius / precision), int((w - ball_radius) / precision) + 1)]
                z_coords = [i * precision for i in range(int(ball_radius / precision), int((h - ball_radius) / precision) + 1)]
                
                # Pruning Step 2: If the box is too small in any dimension to fit even one ball.
                if not x_coords or not y_coords or not z_coords:
                    continue

                # Create a list of candidate center points.
                # Sorting by z, then y, then x ensures a deterministic greedy choice,
                # filling the box from the bottom-front-left corner.
                candidates = []
                for z in z_coords:
                    for y in y_coords:
                        for x in x_coords:
                            candidates.append((x, y, z))
                
                placed_balls_centers = []
                
                while candidates:
                    # Place a ball at the first available candidate position.
                    center = candidates[0]
                    placed_balls_centers.append(center)
                    
                    # Filter the list of candidates, removing any that are too close to the new ball.
                    remaining_candidates = []
                    for cand in candidates[1:]:
                        dist_sq = (center[0] - cand[0])**2 + (center[1] - cand[1])**2 + (center[2] - cand[2])**2
                        if dist_sq >= min_dist_sq:
                            remaining_candidates.append(cand)
                    candidates = remaining_candidates

                num_balls_packed = len(placed_balls_centers)

                # Check if this box is a valid and better solution.
                if num_balls_packed >= min_balls_required:
                    if current_surface_area < min_surface_area:
                        min_surface_area = current_surface_area
                        best_box_dims = (l, w, h)

    # --- Output the result ---
    if best_box_dims:
        l, w, h = best_box_dims
        # The final equation is the set of numbers representing the solution.
        print(f"{l}:{w}:{h}:{min_surface_area}")
    else:
        print("0")

solve_pioneer_packing()
<<<11:11:10:682>>>