import math

def solve_design_problem():
    """
    Solves the container design optimization problem to find the minimum total cost.
    """
    # 1. Define constants from the problem
    r_ball = 2.0  # cm
    d_ball = 4.0  # cm
    energy_per_ball = 30  # MJ
    cost_per_ball = 1000  # USD
    energy_needed = 1000  # MJ
    max_sa = 1000.0  # cm^2
    cost_per_sa = 200  # USD/cm^2
    precision = 0.5  # cm

    # 2. Calculate minimum number of balls
    min_balls = math.ceil(energy_needed / energy_per_ball)

    min_total_cost = float('inf')
    best_config = None

    # 3. Iterate on N, the number of balls, starting from the minimum required
    # Since cost is monotonic with N, the first N that yields a valid solution is the optimum.
    for n_balls in range(min_balls, min_balls + 10):
        ball_cost = n_balls * cost_per_ball

        # If ball_cost alone is more than the best solution found, we can stop.
        if ball_cost > min_total_cost:
            break

        # --- Check Box Container ---
        min_sa_box = float('inf')
        best_box_dims = None
        # To hold n_balls, we need nx*ny*nz >= n_balls
        # Iterate through possible numbers of balls along each axis
        for nx in range(1, n_balls + 1):
            for ny in range(nx, n_balls + 1):
                if nx * ny > n_balls:
                    break
                nz = math.ceil(n_balls / (nx * ny))
                
                l, w, h = nx * d_ball, ny * d_ball, nz * d_ball
                sa = 2 * (l*w + w*h + h*l)

                if sa < min_sa_box:
                    min_sa_box = sa
                    best_box_dims = (l, w, h)

        if min_sa_box <= max_sa:
            container_cost = min_sa_box * cost_per_sa
            total_cost = container_cost + ball_cost
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_config = {
                    "type": "Box", "N": n_balls, "cost": total_cost,
                    "SA": min_sa_box, "dims": best_box_dims,
                    "container_cost": container_cost, "ball_cost": ball_cost
                }

        # --- Check Cylinder Container ---
        min_sa_cyl = float('inf')
        best_cyl_dims = None
        # Packing density for hexagonal circle packing
        packing_density = math.pi / (2 * math.sqrt(3))

        # Iterate through number of layers
        for nz_layers in range(1, n_balls + 1):
            balls_per_layer = math.ceil(n_balls / nz_layers)
            
            # Area of container cross-section = Area of balls / packing_density
            # pi * R^2 = (balls_per_layer * pi * r_ball^2) / packing_density
            r_needed_sq = (balls_per_layer * r_ball**2) / packing_density
            r_needed = math.sqrt(r_needed_sq)
            
            # Quantize radius and height to the required precision
            container_r = math.ceil(r_needed / precision) * precision
            container_h = nz_layers * d_ball

            sa = 2 * math.pi * container_r**2 + 2 * math.pi * container_r * container_h
            
            if sa < min_sa_cyl:
                min_sa_cyl = sa
                best_cyl_dims = (container_r, container_h)

        if min_sa_cyl <= max_sa:
            container_cost = min_sa_cyl * cost_per_sa
            total_cost = container_cost + ball_cost
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_config = {
                    "type": "Cylinder", "N": n_balls, "cost": total_cost,
                    "SA": min_sa_cyl, "dims": best_cyl_dims,
                    "container_cost": container_cost, "ball_cost": ball_cost
                }
        
        # If we found a solution for this N, it's the optimal one. Break the loop.
        if best_config and best_config["N"] == n_balls:
            break

    # 4. Output the result
    if best_config:
        print("Found an optimal container design.")
        print(f"Container Type: {best_config['type']}")
        print(f"Number of Energy Balls: {best_config['N']}")
        if best_config['type'] == 'Cylinder':
            print(f"Dimensions (Radius, Height): {best_config['dims'][0]:.1f} cm, {best_config['dims'][1]:.1f} cm")
        else:
            print(f"Dimensions (L, W, H): {best_config['dims'][0]:.1f} cm, {best_config['dims'][1]:.1f} cm, {best_config['dims'][2]:.1f} cm")
        print(f"Surface Area: {best_config['SA']:.2f} cm^2")
        print("\nCost Calculation:")
        print(f"Total Cost = Container Cost + Ball Cost")
        # Final printout uses rounded values for clarity
        print(f"C = (Surface Area * Cost per cm^2) + (Number of Balls * Cost per Ball)")
        print(f"C = ({best_config['SA']:.2f} * {cost_per_sa}) + ({best_config['N']} * {cost_per_ball})")
        print(f"C = {best_config['container_cost']:.0f} + {best_config['ball_cost']:.0f}")
        print(f"C = {best_config['cost']:.0f}")
        final_cost = round(best_config['cost'])
    else:
        print("Could not find a valid container design with the given constraints.")
        final_cost = 0

    return final_cost

# Execute the function and capture the final answer
final_answer = solve_design_problem()
print(f"\n<<<{final_answer}>>>")
