import math

def solve_design_problem():
    """
    Solves the container design optimization problem.
    Finds the lowest cost container (box or cylinder) to hold a required number of energy balls.
    """

    # 1. Define constants and requirements
    ENERGY_PER_BALL = 30  # MJ
    REQUIRED_ENERGY = 1000  # MJ
    COST_PER_BALL = 1000  # USD
    BALL_RADIUS = 2  # cm
    BALL_DIAMETER = 4  # cm
    MAX_SURFACE_AREA = 1000  # cm^2
    MATERIAL_COST_PER_CM2 = 200  # USD/cm^2
    PRECISION = 0.5  # cm

    # 2. Calculate minimum balls needed and the fixed cost of balls
    n_balls_needed = math.ceil(REQUIRED_ENERGY / ENERGY_PER_BALL)
    cost_of_balls = n_balls_needed * COST_PER_BALL

    min_total_cost = float('inf')
    best_design = None

    # --- Helper function for cylinder packing ---
    memo_balls_per_layer = {}
    def calculate_balls_per_layer(R):
        """Calculates how many balls can be packed in a single circular layer of radius R."""
        if R in memo_balls_per_layer:
            return memo_balls_per_layer[R]
        
        if R < BALL_RADIUS:
            return 0

        # Pack in rows. Start with the central row.
        n_per_layer = math.floor(2 * R / BALL_DIAMETER)
        
        # Add rows offset from the center
        y = BALL_DIAMETER
        while y <= R - BALL_RADIUS: # Center of ball must be inside R-r_ball
            chord_len = 2 * math.sqrt(R**2 - y**2)
            n_in_row = math.floor(chord_len / BALL_DIAMETER)
            if n_in_row > 0:
                n_per_layer += 2 * n_in_row # Add for rows above and below center
            else:
                break
            y += BALL_DIAMETER
        
        memo_balls_per_layer[R] = n_per_layer
        return n_per_layer

    # 3. Search for the optimal Box container
    # Max dimension L for a cube: 6*L^2=1000 -> L~12.9. Max integer step: 12.9/0.5 = 26
    # We need to fit 34 balls, e.g., 4x3x3 grid -> L=16, W=12, H=12.
    # So we need to search up to L=16/0.5=32. Let's give a margin.
    max_dim_steps = int(20 / PRECISION) + 1 
    for l_step in range(1, max_dim_steps):
        L = l_step * PRECISION
        for w_step in range(1, l_step + 1):
            W = w_step * PRECISION
            for h_step in range(1, w_step + 1):
                H = h_step * PRECISION

                surface_area = 2 * (L * W + W * H + H * L)
                if surface_area > MAX_SURFACE_AREA:
                    continue

                balls_fit = math.floor(L / BALL_DIAMETER) * math.floor(W / BALL_DIAMETER) * math.floor(H / BALL_DIAMETER)

                if balls_fit >= n_balls_needed:
                    container_cost = surface_area * MATERIAL_COST_PER_CM2
                    total_cost = container_cost + cost_of_balls
                    if total_cost < min_total_cost:
                        min_total_cost = total_cost
                        best_design = {
                            'type': 'Box',
                            'dims': f"L={L}cm, W={W}cm, H={H}cm",
                            'area': surface_area,
                            'container_cost': container_cost,
                            'total_cost': total_cost,
                            'balls_fit': balls_fit
                        }

    # 4. Search for the optimal Cylinder container
    # Max R for H=R: 4*pi*R^2=1000 -> R~8.9. Max r_step: 8.9/0.5=18
    # Max H for R=2: 2*pi*4 + 2*pi*2*H=1000 -> 4*pi*H=1000-8pi -> H~77. Max h_step: 77/0.5=154
    max_r_steps = int(15 / PRECISION)
    max_h_steps = int(80 / PRECISION)
    for r_step in range(1, max_r_steps):
        R = r_step * PRECISION
        for h_step in range(1, max_h_steps):
            H = h_step * PRECISION

            surface_area = 2 * math.pi * R**2 + 2 * math.pi * R * H
            if surface_area > MAX_SURFACE_AREA:
                continue

            n_layers = math.floor(H / BALL_DIAMETER)
            if n_layers == 0:
                continue
            
            n_per_layer = calculate_balls_per_layer(R)
            balls_fit = n_layers * n_per_layer

            if balls_fit >= n_balls_needed:
                container_cost = surface_area * MATERIAL_COST_PER_CM2
                total_cost = container_cost + cost_of_balls
                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        'type': 'Cylinder',
                        'dims': f"R={R}cm, H={H}cm",
                        'area': surface_area,
                        'container_cost': container_cost,
                        'total_cost': total_cost,
                        'balls_fit': balls_fit
                    }

    # 5. Output the result
    if best_design:
        print(f"Optimal design is a {best_design['type']} with dimensions {best_design['dims']}.")
        print(f"It can hold {best_design['balls_fit']} balls (minimum required is {n_balls_needed}).")
        print(f"Its surface area is {best_design['area']:.2f} cm^2.")
        print("\n--- Cost Breakdown ---")
        print(f"Container Cost = Surface Area * Cost/cm^2")
        print(f"             = {best_design['area']:.2f} cm^2 * {MATERIAL_COST_PER_CM2} USD/cm^2 = {best_design['container_cost']:.2f} USD")
        print(f"Energy Ball Cost = Num Balls * Cost/Ball")
        print(f"                 = {n_balls_needed} * {COST_PER_BALL} USD = {cost_of_balls:.2f} USD")
        print(f"Total Cost (C) = Container Cost + Energy Ball Cost")
        print(f"               = {best_design['container_cost']:.2f} USD + {cost_of_balls:.2f} USD = {best_design['total_cost']:.2f} USD")
        # The final answer C is printed here for clarity, before the special format string.
        print(f"\nFinal Answer C = {best_design['total_cost']:.2f}")

    else:
        print("No solution found that meets all constraints.")
        print("0")

if __name__ == '__main__':
    solve_design_problem()