import math

def solve():
    """
    Solves the container design optimization problem.
    """

    # --- Constants ---
    ENERGY_PER_BALL_MJ = 25
    COST_PER_BALL_USD = 1000
    COST_PER_CM2_USD = 200
    TARGET_ENERGY_MJ = 1000
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    PRECISION_CM = 0.5

    # --- Step 1: Minimum balls needed ---
    min_balls_needed = math.ceil(TARGET_ENERGY_MJ / ENERGY_PER_BALL_MJ)

    # --- Step 2: Box Calculation ---
    def calculate_box_cost():
        min_total_cost_box = float('inf')
        best_config_box = {}

        # Search ball counts from 40 up to a reasonable limit (e.g., 60)
        for n_balls in range(min_balls_needed, min_balls_needed + 20):
            # Find the factors of n_balls that make the most cube-like shape
            best_sa_for_n = float('inf')
            best_factors = None

            for i in range(1, int(n_balls**(1/3)) + 2):
                if n_balls % i == 0:
                    for j in range(i, int((n_balls/i)**0.5) + 2):
                        if (n_balls/i) % j == 0:
                            k = n_balls // (i*j)
                            if i*j*k == n_balls:
                                nx, ny, nz = sorted((i, j, k))
                                
                                L = nx * BALL_DIAMETER_CM
                                W = ny * BALL_DIAMETER_CM
                                H = nz * BALL_DIAMETER_CM
                                
                                surface_area = 2 * (L*W + L*H + W*H)
                                
                                if surface_area < best_sa_for_n:
                                    best_sa_for_n = surface_area
                                    best_factors = (nx, ny, nz)
            
            if best_factors:
                total_cost = (n_balls * COST_PER_BALL_USD) + (best_sa_for_n * COST_PER_CM2_USD)
                
                if total_cost < min_total_cost_box:
                    min_total_cost_box = total_cost
                    best_config_box = {
                        'type': 'Box', 'n_balls': n_balls, 'arrangement': best_factors,
                        'surface_area': best_sa_for_n, 'total_cost': total_cost
                    }
        return best_config_box

    # --- Step 3: Cylinder Calculation ---
    def calculate_cylinder_cost():
        # Packing ratios for fitting k unit-radius circles into a larger circle.
        # k: min_R_factor (where R_container = r_ball * R_factor)
        # Sourced from optimal circle packing data.
        packings = {
            1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0, 8: 3.304,
            9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236, 14: 4.328, 15: 4.521
        }
        
        min_total_cost_cyl = float('inf')
        best_config_cyl = {}

        for k, r_factor in packings.items():
            min_R_math = r_factor * BALL_RADIUS_CM
            # Round up to the nearest multiple of precision
            R_cyl = math.ceil(min_R_math / PRECISION_CM) * PRECISION_CM
            
            n_layers = math.ceil(min_balls_needed / k)
            n_balls = n_layers * k
            H_cyl = n_layers * BALL_DIAMETER_CM
            
            surface_area = 2 * math.pi * R_cyl**2 + 2 * math.pi * R_cyl * H_cyl
            total_cost = (n_balls * COST_PER_BALL_USD) + (surface_area * COST_PER_CM2_USD)
            
            if total_cost < min_total_cost_cyl:
                min_total_cost_cyl = total_cost
                best_config_cyl = {
                    'type': 'Cylinder', 'balls_per_layer': k, 'n_layers': n_layers, 'n_balls': n_balls,
                    'radius': R_cyl, 'height': H_cyl, 'surface_area': surface_area, 'total_cost': total_cost
                }
        return best_config_cyl

    # --- Step 4: Compare and Select Best Design ---
    box_config = calculate_box_cost()
    cyl_config = calculate_cylinder_cost()

    if not box_config and not cyl_config:
        print("Could not find a valid solution.")
        final_cost = 0
    elif not box_config or (cyl_config and cyl_config['total_cost'] < box_config['total_cost']):
        final_config = cyl_config
        # Print final calculation for cylinder
        cost_balls = final_config['n_balls'] * COST_PER_BALL_USD
        cost_material = final_config['surface_area'] * COST_PER_CM2_USD
        print("The optimal design is a Cylinder container.")
        print(f"It holds {final_config['n_balls']} balls, arranged in {final_config['n_layers']} layers with {final_config['balls_per_layer']} balls per layer.")
        print(f"Cylinder Dimensions: Radius = {final_config['radius']} cm, Height = {final_config['height']} cm.")
        print("\nFinal Cost Calculation:")
        print(f"Total Cost = (Number of Balls * Cost per Ball) + (Surface Area * Material Cost per cm^2)")
        print(f"C = ({final_config['n_balls']} * {COST_PER_BALL_USD}) + ((2 * pi * {final_config['radius']}^2 + 2 * pi * {final_config['radius']} * {final_config['height']}) * {COST_PER_CM2_USD})")
        print(f"C = {cost_balls} + ({final_config['surface_area']:.4f} * {COST_PER_CM2_USD})")
        print(f"C = {cost_balls} + {cost_material:.4f}")
        print(f"C = {final_config['total_cost']:.4f}")
        final_cost = final_config['total_cost']
    else:
        final_config = box_config
        # Print final calculation for box
        nx, ny, nz = final_config['arrangement']
        L, W, H = nx * BALL_DIAMETER_CM, ny * BALL_DIAMETER_CM, nz * BALL_DIAMETER_CM
        cost_balls = final_config['n_balls'] * COST_PER_BALL_USD
        cost_material = final_config['surface_area'] * COST_PER_CM2_USD
        print("The optimal design is a Box container.")
        print(f"It holds {final_config['n_balls']} balls, arranged in a {nx}x{ny}x{nz} grid.")
        print(f"Box Dimensions: Length = {L} cm, Width = {W} cm, Height = {H} cm.")
        print("\nFinal Cost Calculation:")
        print(f"Total Cost = (Number of Balls * Cost per Ball) + (Surface Area * Material Cost per cm^2)")
        print(f"C = ({final_config['n_balls']} * {COST_PER_BALL_USD}) + ((2 * ({L}*{W} + {L}*{H} + {W}*{H})) * {COST_PER_CM2_USD})")
        print(f"C = {cost_balls} + ({final_config['surface_area']:.4f} * {COST_PER_CM2_USD})")
        print(f"C = {cost_balls} + {cost_material:.4f}")
        print(f"C = {final_config['total_cost']:.4f}")
        final_cost = final_config['total_cost']
    
    # Final answer in requested format
    return final_cost

# Execute the solution
final_answer = solve()
print(f'<<<{final_answer:.2f}>>>')
