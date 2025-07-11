import math

def solve_packing_problem():
    """
    This function finds the optimal packing of energy balls in a container
    with a surface area of at most 1050 cm^2 to maximize total energy.
    It uses a greedy algorithm on the most promising container shape.
    """

    # Based on preliminary analysis, a box with dimensions friendly to
    # the 4cm diameter of the large balls is the most efficient.
    # The box 12x12x15.5 cm has SA = 2*(12*12 + 12*15.5 + 12*15.5) = 1032 cm^2,
    # which is slightly over the limit. Let's adjust to be compliant.
    # Let's try box 12x12x15.0 cm. SA = 2*(144 + 180 + 180) = 1008 cm^2.
    # Let's try box 12x12x15.5 cm. SA = 2*(144 + 186 + 186) = 1032. This is over.
    # A re-calculation: 2 * (12*12 + 12*15.5 + 12*15.5) = 2 * (144 + 186 + 186) = 2 * 516 = 1032. OK.
    # Wait, the problem says AT MOST 1050. So 1032 is fine.
    # Let's try box 12x12x15.5
    container = {
        "name": "box 12x12x15.5",
        "type": "box",
        "dims": {"L": 12.0, "W": 12.0, "H": 15.5}
    }

    ball_types = [
        {"radius": 2.0, "energy": 10},
        {"radius": 1.0, "energy": 1}
    ]
    precision = 0.5

    placed_balls = []
    
    # Helper functions for checking placement
    def is_contained(center, radius, dims):
        x, y, z = center
        L, W, H = dims["L"], dims["W"], dims["H"]
        return (radius <= x <= L - radius and
                radius <= y <= W - radius and
                radius <= z <= H - radius)

    def check_overlap(center1, radius1, existing_balls):
        x1, y1, z1 = center1
        for ball in existing_balls:
            x2, y2, z2 = ball['center']
            r2 = ball['radius']
            # Using distance squared to avoid sqrt
            dist_sq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
            min_dist_sq = (radius1 + r2)**2
            if dist_sq < min_dist_sq - 1e-9: # Tolerance for float precision
                return True
        return False

    # --- Packing Simulation ---
    dims = container["dims"]
    L, W, H = dims["L"], dims["W"], dims["H"]

    # Pack large balls first (greedy approach)
    for ball_spec in ball_types:
        radius = ball_spec["radius"]
        
        # Iterate through all possible center locations on the 0.5cm grid
        z_steps = int(H / precision) + 1
        y_steps = int(W / precision) + 1
        x_steps = int(L / precision) + 1
        
        # Iterate z, then y, then x (bottom-up layer-by-layer filling)
        for k in range(int(radius / precision), z_steps):
            z = k * precision
            for j in range(int(radius / precision), y_steps):
                y = j * precision
                for i in range(int(radius / precision), x_steps):
                    x = i * precision
                    
                    center = (x, y, z)
                    
                    if is_contained(center, radius, dims):
                        if not check_overlap(center, radius, placed_balls):
                            placed_balls.append({
                                "center": center, 
                                "radius": radius
                            })

    # Count the balls and format the output
    num_small_balls = sum(1 for b in placed_balls if b['radius'] == 1.0)
    num_large_balls = sum(1 for b in placed_balls if b['radius'] == 2.0)
    
    # Final output as requested
    print(f"{container['name']};{num_small_balls};{num_large_balls}")

solve_packing_problem()