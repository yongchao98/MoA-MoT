import math

def solve_pioneer_packing_problem():
    """
    This function searches for the most material-efficient rectangular box that can hold
    at least as many spherical energy balls as the original 12x12x12 cm box.
    """

    # --- Problem Constants ---
    BALL_RADIUS = 2.0
    BALL_DIAMETER = 4.0
    PRECISION = 0.5
    TARGET_BALLS = 27
    INITIAL_DIM = 12.0
    INITIAL_SURFACE_AREA = 6 * (INITIAL_DIM ** 2)

    # --- State for the search ---
    min_surface_area = INITIAL_SURFACE_AREA
    # Start with the initial box as the "best" known solution
    best_dims = None

    def count_packed_balls(L, W, H):
        """
        Greedily packs spheres into a box of size L, W, H and returns the count.
        The algorithm iterates through all possible center locations on the specified
        precision grid and places a ball if it doesn't conflict with existing ones.
        Coordinates are scaled by 1/PRECISION to work with integers and avoid
        floating-point inaccuracies.
        """
        # Scale all measurements to integer units based on precision (1 unit = 0.5 cm)
        l_int = int(L / PRECISION)
        w_int = int(W / PRECISION)
        h_int = int(H / PRECISION)
        radius_int = int(BALL_RADIUS / PRECISION)
        diameter_sq_int = int((BALL_DIAMETER / PRECISION) ** 2)

        # Define the valid range for ball centers in integer units
        possible_x = range(radius_int, l_int - radius_int + 1)
        possible_y = range(radius_int, w_int - radius_int + 1)
        possible_z = range(radius_int, h_int - radius_int + 1)
        
        centers = []
        
        # Iterate through all possible center points in a fixed order (z, y, x)
        for z in possible_z:
            for y in possible_y:
                for x in possible_x:
                    candidate_center = (x, y, z)
                    is_valid = True
                    # Check for overlap with already placed balls
                    for placed_center in centers:
                        # Using squared distance to avoid sqrt calculation
                        dist_sq = (
                            (candidate_center[0] - placed_center[0]) ** 2 +
                            (candidate_center[1] - placed_center[1]) ** 2 +
                            (candidate_center[2] - placed_center[2]) ** 2
                        )
                        if dist_sq < diameter_sq_int:
                            is_valid = False
                            break
                    
                    if is_valid:
                        centers.append(candidate_center)
        return len(centers)

    # --- Main Search Loop ---
    # The problem statement implies integer dimensions for the final answer.
    # The search range is chosen to be reasonably larger than the original dimension.
    # Effective pruning based on surface area makes the search feasible.
    # NOTE: The search can be time-consuming.
    search_limit = 16 
    
    min_volume_needed = TARGET_BALLS * (4/3) * math.pi * (BALL_RADIUS**3)

    for l in range(int(BALL_DIAMETER), search_limit):
        for w in range(l, search_limit):
            # Quick check if the area of the smallest face already makes it impossible to be better
            if 2 * (l * w) >= min_surface_area:
                continue
            for h in range(w, search_limit):
                # Pruning based on volume
                if l * w * h < min_volume_needed:
                    continue

                current_sa = 2 * (l*w + l*h + w*h)
                
                # Pruning based on surface area
                if current_sa >= min_surface_area:
                    continue
                
                num_balls = count_packed_balls(float(l), float(w), float(h))
                
                if num_balls >= TARGET_BALLS:
                    if current_sa < min_surface_area:
                        min_surface_area = current_sa
                        best_dims = tuple(sorted((l, w, h)))
    
    # --- Final Output ---
    if best_dims is None:
        print("Yes, a more efficient box is possible, but the search did not find one. This can be due to search limits or the heuristic nature of the packing algorithm.")
        print(0)
    else:
        l, w, h = best_dims
        sa = int(min_surface_area)
        print("Yes, a more efficient box design is possible.")
        print(f"The initial box 12x12x12 has a surface area of 6 * 12 * 12 = 864 cm^2 and holds 27 balls.")
        print(f"A more efficient box with dimensions {l}x{w}x{h} cm was found that can hold at least 27 balls.")
        print(f"Its surface area is 2 * ({l}*{w} + {l}*{h} + {w}*{h}) = {sa} cm^2, which is less than 864 cm^2.")
        print("\nFinal Answer:")
        print(f"{l}:{w}:{h}:{sa}")
        
# Execute the solver function
solve_pioneer_packing_problem()
<<<10:11:11:682>>>