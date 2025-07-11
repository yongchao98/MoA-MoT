import math

def find_best_design():
    """
    Solves the container design optimization problem to find the minimum total cost.
    """
    # Problem Parameters
    ENERGY_REQUIRED = 1000  # MJ
    ENERGY_PER_BALL = 25  # MJ
    BALL_COST = 1000  # usd
    MATERIAL_COST = 200  # usd per cm^2
    BALL_DIAMETER = 4.0 # cm
    PRECISION = 0.5 # cm

    # Step 1: Calculate minimum number of balls
    min_balls_needed = math.ceil(ENERGY_REQUIRED / ENERGY_PER_BALL)

    min_total_cost = float('inf')
    best_design = {}

    # --- Box Optimization ---
    # Search for an optimal number of balls to pack, up to double the minimum requirement.
    # A slightly larger number might pack into a more efficient, cube-like shape.
    for n_balls_candidate in range(min_balls_needed, min_balls_needed * 2):
        # Find the most "cube-like" factor triple (nx, ny, nz) for the candidate
        # number of balls, as this configuration minimizes surface area for a given volume.
        best_surface_area_term = float('inf')
        optimal_factors = (0, 0, 0)

        for i in range(1, int(n_balls_candidate**(1/3.0)) + 2):
            if n_balls_candidate % i == 0:
                for j in range(i, int((n_balls_candidate/i)**(1/2.0)) + 2):
                    if (n_balls_candidate / i) % j == 0:
                        k = n_balls_candidate // (i * j)
                        if i * j * k == n_balls_candidate: # ensure we have an exact factorization
                            surface_area_term = i*j + j*k + k*i
                            if surface_area_term < best_surface_area_term:
                                best_surface_area_term = surface_area_term
                                optimal_factors = (i, j, k)
        
        if optimal_factors == (0, 0, 0): continue
        
        nx, ny, nz = optimal_factors
        
        # Dimensions are based on ball arrangement. 4cm is a multiple of 0.5cm.
        L = nx * BALL_DIAMETER
        W = ny * BALL_DIAMETER
        H = nz * BALL_DIAMETER

        surface_area = 2 * (L*W + L*H + W*H)
        container_cost = surface_area * MATERIAL_COST
        balls_cost = n_balls_candidate * BALL_COST
        total_cost = balls_cost + container_cost
        
        if total_cost < min_total_cost:
            min_total_cost = total_cost
            best_design = {
                'type': 'Box', 'cost': total_cost, 'balls': n_balls_candidate,
                'dims': (L, W, H), 'surface_area': surface_area,
                'ball_cost': balls_cost, 'container_cost': container_cost
            }
            
    # --- Cylinder Optimization ---
    # We assume balls are packed in stacked rectangular layers.
    # Iterate through a reasonable number of possible stacked layers.
    for num_layers in range(1, min_balls_needed + 2):
        # Calculate balls needed per layer
        balls_per_layer_req = math.ceil(min_balls_needed / num_layers)

        # Pack into a rectangular grid that is as square as possible
        nx_layer = math.ceil(math.sqrt(balls_per_layer_req))
        ny_layer = math.ceil(balls_per_layer_req / nx_layer)
        
        total_balls_packed = nx_layer * ny_layer * num_layers
        if total_balls_packed < min_balls_needed: continue

        # The rectangular grid of balls has these dimensions
        layer_L = nx_layer * BALL_DIAMETER
        layer_W = ny_layer * BALL_DIAMETER
        
        # Find the smallest cylinder radius to fit this grid, respecting precision
        min_radius = math.sqrt(layer_L**2 + layer_W**2) / 2.0
        R = math.ceil(min_radius / PRECISION) * PRECISION # Round up to nearest 0.5cm
        H = num_layers * BALL_DIAMETER
        
        surface_area = (2 * math.pi * R * H) + (2 * math.pi * R**2)
        container_cost = surface_area * MATERIAL_COST
        balls_cost = total_balls_packed * BALL_COST
        total_cost = balls_cost + container_cost
        
        if total_cost < min_total_cost:
            min_total_cost = total_cost
            best_design = {
                'type': 'Cylinder', 'cost': total_cost, 'balls': total_balls_packed,
                'dims': (R, H), 'surface_area': surface_area,
                'ball_cost': balls_cost, 'container_cost': container_cost
            }

    # --- Output Final Result ---
    if not best_design:
        print("No solution found.")
        final_answer = 0
    else:
        print("--- Optimal Design Analysis ---")
        print(f"The container with the lowest total cost is a {best_design['type']}.")
        
        if best_design['type'] == 'Box':
            L, W, H = best_design['dims']
            print(f"The box dimensions are L={L:.1f} cm, W={W:.1f} cm, H={H:.1f} cm.")
        else: # Cylinder
            R, H = best_design['dims']
            print(f"The cylinder dimensions are R={R:.1f} cm, H={H:.1f} cm.")

        print("\n--- Final Cost Calculation ---")
        print(f"The calculation uses {best_design['balls']} balls, each costing ${BALL_COST}.")
        print(f"The container requires {best_design['surface_area']:.2f} cm^2 of material at ${MATERIAL_COST}/cm^2.")
        
        print("\nFinal Equation:")
        print(f"Total Cost = (Number of Balls * Cost per Ball) + (Surface Area * Material Cost per cm^2)")
        print(f"C = ({best_design['balls']} * {BALL_COST}) + ({best_design['surface_area']:.2f} * {MATERIAL_COST})")
        print(f"C = {best_design['ball_cost']} + {best_design['container_cost']:.2f}")
        print(f"C = {best_design['cost']:.2f}")

        final_answer = round(best_design['cost'])

    print(f"\n<<<{final_answer}>>>")

find_best_design()