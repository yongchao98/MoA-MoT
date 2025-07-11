import math

def solve_container_problem():
    """
    Solves the container design optimization problem.
    """
    # 1. Constants and Requirements
    MIN_ENERGY_MJ = 1000
    BALL_ENERGY_MJ = 30
    BALL_COST_USD = 1000
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4

    MATERIAL_COST_PER_CM2 = 200
    MAX_SA_CM2 = 1000
    PRECISION_CM = 0.5

    # Minimum number of balls needed
    min_n_balls = math.ceil(MIN_ENERGY_MJ / BALL_ENERGY_MJ)
    ball_total_cost = min_n_balls * BALL_COST_USD

    min_total_cost = float('inf')
    best_design = None

    # Helper function to generate dimensions
    def drange(start, stop, step):
        r = start
        while r < stop:
            yield r
            r += step

    # 2. Search for the best BOX container
    # Search range for dimensions in cm. Let's assume dimensions won't exceed 40cm.
    # l >= w >= h to avoid duplicate checks. Smallest dimension must hold one ball.
    for h in drange(BALL_DIAMETER_CM, 30, PRECISION_CM):
        for w in drange(h, 30, PRECISION_CM):
            # Optimization: if partial SA is already too large, break
            if 2 * (w*h) > MAX_SA_CM2:
                break
            for l in drange(w, 40, PRECISION_CM):
                sa = 2 * (l * w + w * h + h * l)
                if sa > MAX_SA_CM2:
                    break

                capacity = math.floor(l / BALL_DIAMETER_CM) * \
                           math.floor(w / BALL_DIAMETER_CM) * \
                           math.floor(h / BALL_DIAMETER_CM)

                if capacity >= min_n_balls:
                    total_cost = sa * MATERIAL_COST_PER_CM2 + ball_total_cost
                    if total_cost < min_total_cost:
                        min_total_cost = total_cost
                        best_design = {
                            "type": "box", "l": l, "w": w, "h": h,
                            "sa": sa, "capacity": capacity
                        }

    # 3. Search for the best CYLINDER container
    # Pre-computed data for optimal circle packing (Number of circles -> required R_container/r_circle)
    # Using data from Wikipedia: "Circle packing in a circle"
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.502, 10: 3.693, 11: 3.813, 12: 3.923, 13: 4.029,
        14: 4.236, 15: 4.328, 16: 4.440, 17: 4.521, 18: 4.615, 19: 4.708,
        34: 6.134 # We need at least 34 balls total, but may need fewer per layer
    }
    
    memoized_balls_per_layer = {}
    def get_balls_per_layer(r_container):
        if r_container in memoized_balls_per_layer:
            return memoized_balls_per_layer[r_container]

        ratio = r_container / BALL_RADIUS_CM
        max_n = 0
        # Iterate through the packing data to find the max number of circles that fit
        for n, req_ratio in packing_ratios.items():
            if ratio >= req_ratio:
                max_n = max(max_n, n)
        
        # A more general search for larger N not in the table
        # N ~ 0.9069 * (R/r)^2 for hexagonal packing, can be used as an approximation
        approx_n = math.floor(0.9069 * (ratio ** 2))
        max_n = max(max_n, approx_n)

        memoized_balls_per_layer[r_container] = max_n
        return max_n

    # Search range for dimensions in cm
    for h in drange(BALL_DIAMETER_CM, 40, PRECISION_CM):
        # Find max radius for this height that satisfies SA <= 1000
        # 2*pi*r^2 + 2*pi*r*h - 1000 <= 0. Solve for r.
        c = -MAX_SA_CM2
        b = 2 * math.pi * h
        a = 2 * math.pi
        if (b**2 - 4*a*c) < 0:
            continue
        max_r = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
        
        for r in drange(BALL_RADIUS_CM, max_r + PRECISION_CM, PRECISION_CM):
            sa = 2 * math.pi * r * r + 2 * math.pi * r * h
            if sa > MAX_SA_CM2:
                break
                
            layers = math.floor(h / BALL_DIAMETER_CM)
            if layers == 0:
                continue

            balls_per_layer = get_balls_per_layer(r)
            if balls_per_layer == 0:
                continue
            
            capacity = layers * balls_per_layer

            if capacity >= min_n_balls:
                total_cost = sa * MATERIAL_COST_PER_CM2 + ball_total_cost
                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_design = {
                        "type": "cylinder", "r": r, "h": h,
                        "sa": sa, "capacity": capacity
                    }

    # 4. Output the result
    if best_design:
        sa = best_design['sa']
        container_cost = sa * MATERIAL_COST_PER_CM2
        total_cost = min_total_cost
        
        print("Design Found!")
        if best_design['type'] == 'box':
            print(f"Optimal shape: Box (l={best_design['l']} cm, w={best_design['w']} cm, h={best_design['h']} cm)")
        else:
            print(f"Optimal shape: Cylinder (r={best_design['r']} cm, h={best_design['h']} cm)")
            
        print(f"The container's surface area is {sa:.4f} cm^2.")
        print(f"It can hold {best_design['capacity']} balls (minimum required is {min_n_balls}).")
        
        print("\nCost Calculation:")
        print(f"Container Material Cost = {sa:.4f} cm^2 * {MATERIAL_COST_PER_CM2} usd/cm^2 = {container_cost:.2f} usd")
        print(f"Energy Ball Cost = {min_n_balls} balls * {BALL_COST_USD} usd/ball = {ball_total_cost:.2f} usd")
        print(f"Total Cost = {container_cost:.2f} + {ball_total_cost:.2f} = {total_cost:.2f} usd")
        
        print(f"\nFinal Answer: C = {total_cost:.2f}")
        print(f"<<<{total_cost:.2f}>>>")
    else:
        print("No container design found that meets all criteria.")
        print("C = 0")
        print("<<<0>>>")

solve_container_problem()
