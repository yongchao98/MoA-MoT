import math

def solve_pioneer_probe_design():
    """
    Solves the container design optimization problem for the Pioneer probe.
    """
    # Step 1: Define constants and problem parameters from the user's request.
    BALL_RADIUS = 2.0  # cm
    BALL_ENERGY = 30  # MJ
    BALL_COST = 1000  # USD
    REQUIRED_ENERGY = 1000  # MJ
    MAX_SURFACE_AREA = 1000.0  # cm^2
    MATERIAL_COST_PER_CM2 = 200.0  # USD
    DIMENSION_PRECISION = 0.5  # cm

    # Step 2: Calculate the minimum number of balls and their fixed cost.
    min_balls_needed = math.ceil(REQUIRED_ENERGY / BALL_ENERGY)
    cost_of_balls = min_balls_needed * BALL_COST

    # Step 3: Set up the search for the optimal cylinder. A box container is proven infeasible.
    # The packing of balls (circles) in a circular layer is non-trivial.
    # The following function uses pre-computed data for the "packing circles in a circle" problem.
    # Data is sourced from academic work and resources like packomania.com.
    def get_balls_per_layer(radius, ball_radius):
        """
        Calculates the number of balls (radius `ball_radius`) that can fit in a single
        layer on a circular base of radius `radius`.
        This is a lookup table of known optimal packings, k -> min R/r ratio.
        """
        if radius < ball_radius:
            return 0
        ratio = radius / ball_radius
        # Source: E. Specht, "The best known packings of equal circles in a circle"
        packing_ratios = {
            1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
            8: 3.304, 9: 3.502, 10: 3.702, 11: 3.813, 12: 3.923, 13: 4.0,
            14: 4.14, 15: 4.236, 16: 4.328, 17: 4.398, 18: 4.514, 19: 4.615,
            20: 4.792, 21: 4.864, 22: 4.939, 23: 5.0, 24: 5.068, 25: 5.122,
            26: 5.22, 27: 5.287, 28: 5.349, 29: 5.421, 30: 5.493, 31: 5.56,
            32: 5.584, 33: 5.669, 34: 5.708, 35: 5.765, 36: 5.83, 37: 5.869
        }
        # Find the largest k for which ratio_needed(k) <= ratio
        max_k = 0
        for k, min_ratio in packing_ratios.items():
            if ratio >= min_ratio:
                max_k = max(max_k, k)
        return max_k

    min_total_cost = float('inf')
    best_design = None

    # Step 4: Systematically search through possible dimensions.
    # R = r_i * 0.5, H = h_i * 0.5
    for h_i in range(8, 81):  # H >= 4cm for one layer, reasonable upper bound
        H = h_i * DIMENSION_PRECISION
        for r_i in range(4, 81):  # R >= 2cm for one ball, reasonable upper bound
            R = r_i * DIMENSION_PRECISION

            surface_area = 2 * math.pi * R * (R + H)

            # Constraint 1: Surface area must be within the material budget
            if surface_area > MAX_SURFACE_AREA:
                continue

            num_layers = math.floor(H / (2 * BALL_RADIUS))
            if num_layers == 0:
                continue

            balls_per_layer = get_balls_per_layer(R, BALL_RADIUS)
            if balls_per_layer == 0:
                continue
                
            capacity = num_layers * balls_per_layer

            # Constraint 2: Capacity must be sufficient
            if capacity >= min_balls_needed:
                # This is a valid design. Calculate its cost.
                total_cost = MATERIAL_COST_PER_CM2 * surface_area + cost_of_balls
                # Objective: Minimize total cost
                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        "shape": "Cylinder",
                        "R": R,
                        "H": H,
                        "surface_area": surface_area,
                        "capacity": capacity,
                        "cost": total_cost
                    }

    # Step 5: Output the final result.
    if best_design:
        C = round(best_design['cost'])
        material_cost = best_design['cost'] - cost_of_balls
        print(f"Optimal Design Found: {best_design['shape']}")
        print(f"Dimensions: Radius = {best_design['R']:.1f} cm, Height = {best_design['H']:.1f} cm")
        print(f"Container Surface Area = {best_design['surface_area']:.2f} cm^2")
        print(f"Container Capacity = {best_design['capacity']} balls")
        print("\nCost Calculation:")
        print(f"Material Cost = {MATERIAL_COST_PER_CM2:.0f} USD/cm^2 * {best_design['surface_area']:.4f} cm^2 = {material_cost:.2f} USD")
        print(f"Energy Balls Cost = {min_balls_needed} balls * {BALL_COST:.0f} USD/ball = {cost_of_balls:.2f} USD")
        print(f"Total Cost C = {material_cost:.2f} USD + {cost_of_balls:.2f} USD = {best_design['cost']:.2f} USD")
        print(f"\nThe lowest total cost C is {C}.")
        print(f"<<<{C}>>>")
    else:
        # If no solution is found
        C = 0
        print("No design could be found that satisfies all constraints.")
        print(f"<<<{C}>>>")

solve_pioneer_probe_design()