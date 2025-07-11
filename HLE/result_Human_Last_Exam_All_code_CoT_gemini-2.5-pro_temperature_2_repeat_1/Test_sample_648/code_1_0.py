import math

def solve_pioneer_probe_design():
    """
    Finds the optimal container design by searching through possible dimensions
    for both box and cylinder shapes, while respecting all constraints.
    """

    # --- Step 1: Define constants and initialize tracking variables ---
    BALL_DIAMETER = 4.0
    MIN_BALLS_REQUIRED = 34
    MAX_SURFACE_AREA = 1000.0
    COST_PER_CM2 = 200.0
    COST_PER_BALL = 1000.0
    PRECISION = 0.5

    min_total_cost = float('inf')
    best_design_details = None

    # --- Helper function for calculating balls in a cylinder layer (grid packing) ---
    def get_balls_per_layer(radius, ball_diameter):
        ball_radius = ball_diameter / 2.0
        # The centers of balls must fit within a circle of radius R - r_ball
        effective_radius = radius - ball_radius
        if effective_radius < 0:
            return 0
        
        count = 0
        # We check for ball centers placed on a square grid of size ball_diameter
        y_step_max = math.floor(effective_radius / ball_diameter)
        
        # Count balls in the central row (y=0)
        x_max_at_center = effective_radius
        x_step_max = math.floor(x_max_at_center / ball_diameter)
        count += (2 * x_step_max + 1)
        
        # Count balls in rows above (and symmetrically, below) the center
        for j in range(1, y_step_max + 1):
            y = j * ball_diameter
            if y > effective_radius:
                continue
            x_max_at_y = math.sqrt(effective_radius**2 - y**2)
            x_step_max_at_y = math.floor(x_max_at_y / ball_diameter)
            # Add balls for this row (multiplied by 2 for positive and negative y)
            count += 2 * (2 * x_step_max_at_y + 1)
            
        return count

    # --- Step 2: Search for the optimal Box container ---
    min_dim_steps = int(BALL_DIAMETER / PRECISION)
    max_l_steps = int(60.5 / PRECISION) + 2  # Max length if W,H are minimal

    for l_steps in range(min_dim_steps, max_l_steps):
        L = l_steps * PRECISION
        for w_steps in range(min_dim_steps, l_steps + 1):
            W = w_steps * PRECISION
            for h_steps in range(min_dim_steps, w_steps + 1):
                H = h_steps * PRECISION
                
                surface_area = 2 * (L * W + W * H + H * L)
                if surface_area > MAX_SURFACE_AREA:
                    continue

                balls_packed = math.floor(L / BALL_DIAMETER) * math.floor(W / BALL_DIAMETER) * math.floor(H / BALL_DIAMETER)
                if balls_packed < MIN_BALLS_REQUIRED:
                    continue
                
                total_cost = surface_area * COST_PER_CM2 + balls_packed * COST_PER_BALL
                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design_details = {
                        "sa": surface_area, "balls": balls_packed, "cost": total_cost
                    }

    # --- Step 3: Search for the optimal Cylinder container ---
    min_r_steps = int((BALL_DIAMETER / 2.0) / PRECISION)
    min_h_steps = int(BALL_DIAMETER / PRECISION)
    max_r_steps = int(12.7 / PRECISION) + 1  # Max radius is approx 12.6
    max_h_steps = int(78.0 / PRECISION) + 1  # Max height is approx 77.5
    
    for h_steps in range(min_h_steps, max_h_steps):
        H = h_steps * PRECISION
        for r_steps in range(min_r_steps, max_r_steps):
            R = r_steps * PRECISION
            
            surface_area = 2 * math.pi * R * H + 2 * math.pi * R**2
            if surface_area > MAX_SURFACE_AREA:
                continue

            balls_per_layer = get_balls_per_layer(R, BALL_DIAMETER)
            if balls_per_layer == 0: continue
            
            num_layers = math.floor(H / BALL_DIAMETER)
            if num_layers == 0: continue
                
            balls_packed = balls_per_layer * num_layers
            if balls_packed < MIN_BALLS_REQUIRED:
                continue

            total_cost = surface_area * COST_PER_CM2 + balls_packed * COST_PER_BALL
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design_details = {
                    "sa": surface_area, "balls": balls_packed, "cost": total_cost
                }

    # --- Step 4: Output the result ---
    if best_design_details is None:
        final_answer = 0
        print(final_answer)
    else:
        sa = best_design_details["sa"]
        balls = best_design_details["balls"]
        cost = best_design_details["cost"]
        final_answer = round(cost)
        print(f"Total Cost = {sa:.2f} * {COST_PER_CM2:.0f} + {balls} * {COST_PER_BALL:.0f} = {final_answer}")
    
    return final_answer

# Execute the solution function and format the final answer.
final_result = solve_pioneer_probe_design()
print(f"<<<{final_result}>>>")