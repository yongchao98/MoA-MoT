import math

def solve_pioneer_probe_problem():
    """
    Solves the Pioneer probe container design problem.
    """
    # --- Problem Constants ---
    REQUIRED_ENERGY_MJ = 1000
    BALL_ENERGY_MJ = 30
    BALL_RADIUS_CM = 2.0
    BALL_DIAMETER_CM = 4.0
    MAX_SURFACE_AREA_CM2 = 1000.0
    MATERIAL_COST_PER_CM2 = 200.0
    BALL_COST = 1000.0
    PRECISION_CM = 0.5

    MIN_BALLS = math.ceil(REQUIRED_ENERGY_MJ / BALL_ENERGY_MJ)

    # --- Box Analysis ---
    # As derived in the plan, no box container is feasible because the minimum
    # surface area for any configuration holding >= 34 balls exceeds 1000 cm^2.
    # The most efficient packing (4x3x3=36 balls) requires a SA of 1056 cm^2.

    # --- Cylinder Analysis ---
    
    # Data for optimal packing of N circles in a larger circle.
    # Format: {N: R_container / r_ball}
    # Source: Packomania and other mathematical sources on circle packing.
    PACKING_RATIOS = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.532, 10: 3.639, 11: 3.813, 12: 3.968, 13: 4.0,
        14: 4.161, 15: 4.236, 16: 4.236, 17: 4.236, 18: 4.236, 19: 4.236,
        20: 4.387, 21: 4.470, 22: 4.470
    }

    def get_max_balls_per_layer(radius_cm):
        """Find max number of balls per layer for a given cylinder radius."""
        ratio = radius_cm / BALL_RADIUS_CM
        max_n = 0
        for n, r_ratio in PACKING_RATIOS.items():
            if ratio >= r_ratio:
                max_n = max(max_n, n)
        return max_n

    min_total_cost = float('inf')
    best_config = None

    # Search space for radius R in cm, using the required precision.
    # Smallest R must be > ball radius. Max R is limited by SA.
    # If H=2R (optimal SA/V), 6*pi*R^2 <= 1000 -> R <= 7.28. We search a bit beyond.
    min_r = math.ceil(BALL_RADIUS_CM / PRECISION_CM) * PRECISION_CM
    max_r = 10.0
    
    r_steps = int((max_r - min_r) / PRECISION_CM) + 1
    for i in range(r_steps):
        r_cm = min_r + i * PRECISION_CM
        
        n_per_layer = get_max_balls_per_layer(r_cm)
        if n_per_layer == 0:
            continue

        num_layers = math.ceil(MIN_BALLS / n_per_layer)
        h_cm = num_layers * BALL_DIAMETER_CM

        surface_area = 2 * math.pi * r_cm * (r_cm + h_cm)

        if surface_area <= MAX_SURFACE_AREA_CM2:
            num_balls_fit = n_per_layer * num_layers
            container_cost = surface_area * MATERIAL_COST_PER_CM2
            balls_cost = num_balls_fit * BALL_COST
            total_cost = container_cost + balls_cost

            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_config = {
                    "type": "Cylinder",
                    "R_cm": r_cm,
                    "H_cm": h_cm,
                    "SA_cm2": surface_area,
                    "n_per_layer": n_per_layer,
                    "num_layers": num_layers,
                    "total_balls": num_balls_fit,
                    "cost": total_cost,
                    "container_cost": container_cost,
                    "balls_cost": balls_cost
                }

    if best_config:
        print("--- Optimal Design Found: Cylinder ---")
        print(f"Radius (R): {best_config['R_cm']:.1f} cm")
        print(f"Height (H): {best_config['H_cm']:.1f} cm")
        print("\n--- Cost Calculation ---")
        print(f"Surface Area = 2 * {math.pi:.5f} * {best_config['R_cm']:.1f} * ({best_config['R_cm']:.1f} + {best_config['H_cm']:.1f}) = {best_config['SA_cm2']:.2f} cm^2")
        print(f"Container Material Cost = {best_config['SA_cm2']:.2f} cm^2 * {MATERIAL_COST_PER_CM2:.0f} usd/cm^2 = {best_config['container_cost']:.2f} usd")
        print("\n--- Energy Ball Packing ---")
        print(f"Balls per layer = {best_config['n_per_layer']}")
        print(f"Number of layers = {best_config['num_layers']}")
        print(f"Total balls packed = {best_config['n_per_layer']} * {best_config['num_layers']} = {best_config['total_balls']}")
        print(f"Energy Ball Cost = {best_config['total_balls']} balls * {BALL_COST:.0f} usd/ball = {best_config['balls_cost']:.2f} usd")
        print("\n--- Final Total Cost ---")
        print(f"Total Cost = {best_config['container_cost']:.2f} usd + {best_config['balls_cost']:.2f} usd = {best_config['cost']:.2f} usd")
        
        global C
        C = best_config['cost']
    else:
        print("No feasible design found that meets all constraints.")
        C = 0

solve_pioneer_probe_problem()
print(f'<<<{C:.2f}>>>')