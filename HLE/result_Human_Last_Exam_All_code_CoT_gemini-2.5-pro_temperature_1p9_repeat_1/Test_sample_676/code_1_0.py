import math

def solve_container_problem():
    """
    Calculates the minimum cost to package energy balls in either a box or cylinder.
    """

    # 1. Define constants from the problem description
    BALL_RADIUS = 2.0
    BALL_DIAMETER = 4.0
    BALL_ENERGY = 25.0
    BALL_COST = 1000.0
    TOTAL_ENERGY_REQ = 1000.0
    MATERIAL_COST_PER_CM2 = 200.0
    PRECISION = 0.5
    PI = math.pi

    # Calculate the minimum number of balls required
    min_balls_needed = math.ceil(TOTAL_ENERGY_REQ / BALL_ENERGY)

    # 2. Function to find the best box container
    def find_best_box():
        min_total_cost = float('inf')
        best_config = {}
        limit = int(min_balls_needed) + 1  # A reasonable search limit for one dimension

        for nx in range(1, limit):
            for ny in range(nx, limit):
                # Calculate the required number of layers (nz)
                nz = math.ceil(min_balls_needed / (nx * ny))
                
                # Determine configuration properties
                actual_balls = nx * ny * nz
                l, w, h = nx * BALL_DIAMETER, ny * BALL_DIAMETER, nz * BALL_DIAMETER
                area = 2 * (l * w + l * h + w * h)
                
                # Calculate costs
                ball_cost = actual_balls * BALL_COST
                material_cost = area * MATERIAL_COST_PER_CM2
                total_cost = ball_cost + material_cost

                if total_cost < min_total_cost:
                    min_total_cost = total_cost
                    best_config = {
                        'type': 'Box',
                        'total_cost': total_cost,
                        'ball_cost': ball_cost,
                        'material_cost': material_cost,
                        'area': area,
                        'balls': actual_balls,
                        'ball_config': (nx, ny, nz),
                        'dims': (l, w, h)
                    }
        return best_config

    # 3. Function to find the best cylinder container
    def find_best_cylinder():
        # R_k/r ratios for packing k circles in a circle. Data from packomania.com
        pack_ratios = {
            1: 1.0, 2: 2.0, 3: 2.1547, 4: 2.4142, 5: 2.7013, 6: 3.0, 7: 3.0, 8: 3.3047, 
            9: 3.6131, 10: 3.8130, 11: 3.9238, 12: 4.0296, 13: 4.2360, 14: 4.3283, 15: 4.5212, 
            16: 4.6154, 17: 4.7920, 18: 4.8637, 19: 4.8637, 20: 5.1223, 21: 5.2109, 22: 5.3228,
            23: 5.4309, 24: 5.5451, 25: 5.6568, 26: 5.7587, 27: 5.8679, 28: 5.9754, 29: 6.0772,
            30: 6.1037, 31: 6.2731, 32: 6.3769, 33: 6.4719, 34: 6.5658, 35: 6.6669, 36: 6.7479,
            37: 6.7479, 38: 6.9388, 39: 7.0379, 40: 7.1352
        }
        
        min_total_cost = float('inf')
        best_config = {}

        for k in range(1, int(min_balls_needed) + 1):
            if k not in pack_ratios: continue

            # Calculate cylinder radius, enforcing precision
            min_radius = pack_ratios[k] * BALL_RADIUS
            radius = math.ceil(min_radius / PRECISION) * PRECISION
            
            # Calculate cylinder height
            num_layers = math.ceil(min_balls_needed / k)
            height = num_layers * BALL_DIAMETER
            
            # Determine configuration properties
            actual_balls = k * num_layers
            area = 2 * PI * radius**2 + 2 * PI * radius * height
            
            # Calculate costs
            ball_cost = actual_balls * BALL_COST
            material_cost = area * MATERIAL_COST_PER_CM2
            total_cost = ball_cost + material_cost

            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_config = {
                    'type': 'Cylinder',
                    'total_cost': total_cost,
                    'ball_cost': ball_cost,
                    'material_cost': material_cost,
                    'area': area,
                    'balls': actual_balls,
                    'ball_config': (k, int(num_layers)),  # k balls per layer, for num_layers
                    'dims': (radius, height)
                }
        return best_config

    # 4. Compare designs and print the results
    box_result = find_best_box()
    cylinder_result = find_best_cylinder()

    if box_result['total_cost'] < cylinder_result['total_cost']:
        final_design = box_result
    else:
        final_design = cylinder_result

    print("--- Optimal Container Design ---")
    print(f"The best container is a {final_design['type']}.")
    
    if final_design['type'] == 'Box':
        nx, ny, nz = final_design['ball_config']
        l, w, h = final_design['dims']
        print(f"It is designed to hold a grid of {nx}x{ny}x{nz} energy balls.")
        print(f"Dimensions: Length = {l:.1f} cm, Width = {w:.1f} cm, Height = {h:.1f} cm.")
    else: # Cylinder
        k, num_layers = final_design['ball_config']
        r, h = final_design['dims']
        print(f"It is designed to hold {num_layers} layers of {k} balls each.")
        print(f"Dimensions: Radius = {r:.1f} cm, Height = {h:.1f} cm.")
    
    print(f"\nSurface Area: {final_design['area']:.2f} cm^2")
    
    # Print the breakdown for the final equation as requested
    print("\n--- Cost Breakdown ---")
    print(f"Total Balls: {int(final_design['balls'])}")
    C_balls = final_design['ball_cost']
    C_material = final_design['material_cost']
    C_total = final_design['total_cost']
    
    print(f"Cost of Energy Balls: {int(C_balls)} = {int(final_design['balls'])} balls * ${BALL_COST:.0f}/ball")
    print(f"Cost of Material: {round(C_material)} = {final_design['area']:.2f} cm^2 * ${MATERIAL_COST_PER_CM2:.0f}/cm^2")
    print("\nFinal Equation:")
    print(f"Total Cost (C) = Ball Cost + Material Cost")
    print(f"{round(C_total)} = {int(C_balls)} + {round(C_material)}")
    
    return round(C_total)

# Run the solver and get the final answer C
C = solve_container_problem()