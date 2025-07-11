import math

def solve_container_problem():
    """
    Solves the container design optimization problem to find the minimum cost.
    """
    # Constants
    BALL_RADIUS = 2.0
    BALL_DIAMETER = BALL_RADIUS * 2
    ENERGY_PER_BALL = 30  # MJ
    COST_PER_BALL = 1000  # USD
    REQUIRED_ENERGY = 1000  # MJ
    MAX_SURFACE_AREA = 1000.0  # cm^2
    COST_PER_SA = 200  # USD per cm^2
    PRECISION_STEP = 0.5  # cm

    # Step 1: Calculate minimum number of balls
    min_balls_needed = math.ceil(REQUIRED_ENERGY / ENERGY_PER_BALL)

    min_total_cost = float('inf')
    best_design = None

    # --- Box Container Search ---
    # Search range for dimension steps (l, w, h are multiples of PRECISION_STEP)
    # Max dimension if others are minimal (4cm): 2(4L + 4L + 16) <= 1000 -> 16L + 32 <= 1000 -> 16L <= 968 -> L <= 60.5
    # So max integer step is around 60.5 / 0.5 = 121. We search up to 125.
    # Min dimension must be >= 4.0 cm, so min step is 4.0 / 0.5 = 8.
    for l_step in range(8, 126):
        for w_step in range(8, l_step + 1):
            for h_step in range(8, w_step + 1):
                L = l_step * PRECISION_STEP
                W = w_step * PRECISION_STEP
                H = h_step * PRECISION_STEP

                # Check SA constraint
                surface_area = 2 * (L * W + L * H + W * H)
                if surface_area > MAX_SURFACE_AREA:
                    # Since H is the smallest dim, further increases in H will also fail
                    if h_step == 8:
                        break # Optimization for inner loop
                    continue

                # Check ball count constraint
                num_balls = math.floor(L / BALL_DIAMETER) * \
                            math.floor(W / BALL_DIAMETER) * \
                            math.floor(H / BALL_DIAMETER)
                if num_balls < min_balls_needed:
                    continue

                # Calculate cost
                material_cost = surface_area * COST_PER_SA
                ball_cost = num_balls * COST_PER_BALL
                total_cost = material_cost + ball_cost

                # Update best design if this one is cheaper
                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        "type": "Box",
                        "dims": {"L": L, "W": W, "H": H},
                        "SA": surface_area,
                        "num_balls": num_balls,
                        "material_cost": material_cost,
                        "ball_cost": ball_cost,
                        "total_cost": total_cost,
                    }

    # --- Cylinder Container Search ---
    # Search range: if R=3, H can be up to 50. if R=8.9, H=8.9.
    # r_cyl_max ~ 8.9/0.5=18. h_cyl_max ~ 50/0.5=100.
    # Search r_cyl up to 40, h_cyl up to 100. Min r_cyl = ceil(2*sqrt(2)/0.5)=6. Min h_cyl=8.
    for r_cyl_step in range(6, 41):
        for h_cyl_step in range(8, 101):
            R = r_cyl_step * PRECISION_STEP
            H = h_cyl_step * PRECISION_STEP

            # Check SA constraint
            surface_area = 2 * math.pi * R**2 + 2 * math.pi * R * H
            if surface_area > MAX_SURFACE_AREA:
                continue

            # Check ball count constraint (inscribed square prism method)
            num_layers = math.floor(H / BALL_DIAMETER)
            side_of_square = math.sqrt(2) * R
            balls_per_side = math.floor(side_of_square / BALL_DIAMETER)
            balls_per_layer = balls_per_side**2
            num_balls = num_layers * balls_per_layer

            if num_balls < min_balls_needed:
                continue

            # Calculate cost
            material_cost = surface_area * COST_PER_SA
            ball_cost = num_balls * COST_PER_BALL
            total_cost = material_cost + ball_cost

            # Update best design
            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design = {
                    "type": "Cylinder",
                    "dims": {"R": R, "H": H},
                    "SA": surface_area,
                    "num_balls": num_balls,
                    "material_cost": material_cost,
                    "ball_cost": ball_cost,
                    "total_cost": total_cost,
                }


    # Final Output
    if best_design:
        print("Found an optimal design:")
        print(f"Container Type: {best_design['type']}")
        print(f"Dimensions (cm): {best_design['dims']}")
        print(f"Surface Area: {best_design['SA']:.2f} cm^2")
        print(f"Number of Balls: {best_design['num_balls']}")
        print("\nCost Calculation:")
        print(f"Material Cost = {best_design['SA']:.2f} cm^2 * ${COST_PER_SA}/cm^2 = ${best_design['material_cost']:.2f}")
        print(f"Energy Ball Cost = {best_design['num_balls']} balls * ${COST_PER_BALL}/ball = ${best_design['ball_cost']:.2f}")
        print(f"Total Cost = ${best_design['material_cost']:.2f} + ${best_design['ball_cost']:.2f} = ${best_design['total_cost']:.2f}")
        print(f"<<<{best_design['total_cost']:.0f}>>>")
    else:
        print("No design could be found that satisfies all the constraints.")
        print("<<<0>>>")

solve_container_problem()
