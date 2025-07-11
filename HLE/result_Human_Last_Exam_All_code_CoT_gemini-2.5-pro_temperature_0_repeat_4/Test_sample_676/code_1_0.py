import math

def solve_pioneer_container_problem():
    """
    Solves the container optimization problem for the Pioneer probe.
    It calculates the minimum cost for both box and cylinder containers
    and outputs the details of the best design.
    """

    # --- Constants ---
    ENERGY_PER_BALL = 25  # MJ
    COST_PER_BALL = 1000  # USD
    REQUIRED_ENERGY = 1000  # MJ
    BALL_RADIUS = 2  # cm
    BALL_DIAMETER = 4  # cm
    MATERIAL_COST_PER_CM2 = 200  # USD
    PRECISION = 0.5  # cm

    # --- Step 1: Calculate minimum number of balls ---
    num_balls_needed = math.ceil(REQUIRED_ENERGY / ENERGY_PER_BALL)

    # --- Step 2: Find the optimal Box container ---
    def calculate_box_cost(min_balls):
        min_cost = float('inf')
        best_details = {}
        # Search up to a slightly larger number of balls to find optimal arrangements
        limit = min_balls + 15 
        for nx in range(1, limit):
            for ny in range(1, limit):
                if nx * ny > limit:
                    break
                nz = math.ceil(min_balls / (nx * ny))
                num_balls_actual = nx * ny * nz
                
                L, W, H = nx * BALL_DIAMETER, ny * BALL_DIAMETER, nz * BALL_DIAMETER
                area = 2 * (L*W + L*H + W*H)
                cost = num_balls_actual * COST_PER_BALL + area * MATERIAL_COST_PER_CM2
                
                if cost < min_cost:
                    min_cost = cost
                    best_details = {
                        'cost': cost, 'balls': num_balls_actual, 'area': area,
                        'dims': (L, W, H), 'config': (nx, ny, nz)
                    }
        return min_cost, best_details

    # --- Step 3: Find the optimal Cylinder container ---
    def calculate_cylinder_cost(min_balls):
        min_cost = float('inf')
        best_details = {}
        
        # Packing data for k unit circles in a circle of radius R_p
        # k: R_p (radius of circle containing k unit-radius circles)
        packing_radii = {
            1: 0, 2: 1.0, 3: 1.1547, 4: 1.4142, 5: 1.641, 6: 2.0, 7: 2.0,
            8: 2.3094, 9: 2.4142, 10: 2.7155, 11: 2.8284, 12: 2.97, 13: 3.0, 14: 3.0
        }
        def get_packing_radius(k):
            if k in packing_radii:
                return packing_radii[k]
            # For larger k, use density approximation
            # Area_k_spheres / Area_container ~ 0.9069 (hexagonal packing density)
            # k*pi*r^2 / (pi*R^2) ~ 0.9069 => R/r = R_p ~ sqrt(k/0.9069)
            return math.sqrt(k / 0.9069) if k > 0 else 0

        limit = min_balls + 15
        for n_layers in range(1, limit):
            k_actual = math.ceil(min_balls / n_layers)
            num_balls_actual = n_layers * k_actual
            
            H = n_layers * BALL_DIAMETER
            R_p = get_packing_radius(k_actual)
            R_calc = BALL_RADIUS * (1 + R_p)
            R = math.ceil(R_calc / PRECISION) * PRECISION
            
            area = 2 * math.pi * R**2 + 2 * math.pi * R * H
            cost = num_balls_actual * COST_PER_BALL + area * MATERIAL_COST_PER_CM2
            
            if cost < min_cost:
                min_cost = cost
                best_details = {
                    'cost': cost, 'balls': num_balls_actual, 'area': area,
                    'dims': (R, H), 'config': (n_layers, k_actual)
                }
        return min_cost, best_details

    # --- Step 4: Compare designs and print the result ---
    cost_box, details_box = calculate_box_cost(num_balls_needed)
    cost_cyl, details_cyl = calculate_cylinder_cost(num_balls_needed)

    if cost_box < cost_cyl:
        min_cost = cost_box
        best_design = "Box"
        details = details_box
    else:
        min_cost = cost_cyl
        best_design = "Cylinder"
        details = details_cyl

    print(f"The optimal design is a {best_design} container.")
    
    if best_design == "Box":
        nx, ny, nz = details['config']
        L, W, H = details['dims']
        print(f"It is designed to hold a {nx}x{ny}x{nz} grid of balls.")
        print(f"Dimensions: Length={L}cm, Width={W}cm, Height={H}cm.")
    else: # Cylinder
        n_l, k = details['config']
        R, H = details['dims']
        print(f"It is designed to hold {n_l} layers of {k} balls each.")
        print(f"Dimensions: Radius={R}cm, Height={H}cm.")

    print(f"Total balls packed: {details['balls']}")
    print(f"Container surface area: {details['area']:.4f} cm^2")
    print("\n--- Final Cost Calculation ---")
    
    # Output the final equation with all numbers
    cost_of_balls = details['balls'] * COST_PER_BALL
    cost_of_material = details['area'] * MATERIAL_COST_PER_CM2
    
    print(f"Cost of Balls = {details['balls']} balls * ${COST_PER_BALL}/ball = ${cost_of_balls}")
    print(f"Cost of Material = {details['area']:.4f} cm^2 * ${MATERIAL_COST_PER_CM2}/cm^2 = ${cost_of_material:.4f}")
    print("\nTotal Cost C = Cost of Balls + Cost of Material")
    print(f"C = {cost_of_balls} + {cost_of_material:.4f} = ${min_cost:.4f}")

    final_answer = math.ceil(min_cost)
    print(f"\nThe lowest total cost C, rounded up to the nearest dollar, is ${final_answer}.")
    print(f"<<<{final_answer}>>>")

# Execute the solution
solve_pioneer_container_problem()