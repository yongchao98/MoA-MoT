import math

def solve_container_problem():
    """
    Solves the energy ball container design problem by finding the minimum cost
    for both box and cylinder designs and selecting the best one.
    """
    # --- Problem Parameters ---
    BALL_RADIUS_CM = 2.0
    BALL_DIAMETER_CM = 4.0
    ENERGY_PER_BALL_MJ = 25
    COST_PER_BALL_USD = 1000
    REQUIRED_ENERGY_MJ = 1000
    COST_PER_CM2_USD = 200
    PRECISION_CM = 0.5
    
    MIN_BALLS_REQUIRED = math.ceil(REQUIRED_ENERGY_MJ / ENERGY_PER_BALL_MJ)

    # --- Find Best Box Design ---
    def find_best_box():
        best_box = {'cost': float('inf')}
        
        # Search for number of balls from 40 up to 60.
        # It's unlikely that using many more than 40 balls will be cheaper.
        for num_balls in range(MIN_BALLS_REQUIRED, MIN_BALLS_REQUIRED + 20):
            # Find all integer factor triples (n_l, n_w, n_h) for num_balls
            for nl in range(1, int(num_balls** (1/3)) + 2):
                if num_balls % nl == 0:
                    for nw in range(nl, int((num_balls/nl)**0.5) + 2):
                        if (num_balls / nl) % nw == 0:
                            nh = num_balls // (nl * nw)
                            if nl * nw * nh == num_balls:
                                # Dimensions based on simple cubic packing
                                l = nl * BALL_DIAMETER_CM
                                w = nw * BALL_DIAMETER_CM
                                h = nh * BALL_DIAMETER_CM
                                
                                # Surface Area
                                surface_area = 2 * (l*w + l*h + w*h)
                                
                                # Cost
                                cost_balls = num_balls * COST_PER_BALL_USD
                                cost_container = surface_area * COST_PER_CM2_USD
                                total_cost = cost_balls + cost_container
                                
                                if total_cost < best_box['cost']:
                                    best_box = {
                                        'cost': total_cost,
                                        'type': 'Box',
                                        'balls': num_balls,
                                        'layout': (nl, nw, nh),
                                        'dims': (l, w, h),
                                        'area': surface_area,
                                        'cost_balls': cost_balls,
                                        'cost_container': cost_container
                                    }
        return best_box

    # --- Find Best Cylinder Design ---
    def find_best_cylinder():
        best_cylinder = {'cost': float('inf')}
        
        # Packing ratio C(n) = R_container / r_ball for n circles in a circle.
        # Data from packomania.com for optimal packings.
        C_n_lookup = {
            1: 1.0, 2: 2.0, 3: 2.154, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
            8: 3.304, 9: 3.534, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236,
            14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.863, 19: 5.0,
            20: 5.122, 21: 5.248, 22: 5.336, 23: 5.485, 24: 5.568, 25: 5.679,
            26: 5.760, 27: 5.864, 28: 5.965, 29: 6.090, 30: 6.155, 31: 6.273,
            32: 6.347, 33: 6.425, 34: 6.502, 35: 6.579, 36: 6.657, 37: 6.747,
            38: 6.811, 39: 6.903, 40: 6.941
        }
        
        # Iterate through number of balls per layer
        for n_per_layer in range(1, MIN_BALLS_REQUIRED + 1):
            if n_per_layer not in C_n_lookup:
                continue

            num_layers = math.ceil(MIN_BALLS_REQUIRED / n_per_layer)
            num_balls = n_per_layer * num_layers
            
            # Calculate Radius R, applying precision rule
            c_n = C_n_lookup[n_per_layer]
            radius_raw = c_n * BALL_RADIUS_CM
            radius = math.ceil(radius_raw / PRECISION_CM) * PRECISION_CM
            
            # Calculate Height H using dense packing, applying precision rule
            if num_layers == 1:
                height_raw = BALL_DIAMETER_CM
            else:
                # Dense stacking (HCP/FCC), layers nestle into each other
                height_raw = BALL_DIAMETER_CM + (num_layers - 1) * BALL_DIAMETER_CM * math.sqrt(2/3)
            height = math.ceil(height_raw / PRECISION_CM) * PRECISION_CM

            # Surface Area
            surface_area = (2 * math.pi * radius**2) + (2 * math.pi * radius * height)
            
            # Cost
            cost_balls = num_balls * COST_PER_BALL_USD
            cost_container = surface_area * COST_PER_CM2_USD
            total_cost = cost_balls + cost_container

            if total_cost < best_cylinder['cost']:
                best_cylinder = {
                    'cost': total_cost,
                    'type': 'Cylinder',
                    'balls': num_balls,
                    'layout': f"{num_layers} layers of {n_per_layer} balls",
                    'dims': (radius, height),
                    'area': surface_area,
                    'cost_balls': cost_balls,
                    'cost_container': cost_container
                }
        return best_cylinder

    # --- Determine the optimal design ---
    best_box = find_best_box()
    best_cylinder = find_best_cylinder()

    if best_box['cost'] < best_cylinder['cost']:
        best_design = best_box
    else:
        best_design = best_cylinder
        
    final_cost = round(best_design['cost'])

    # --- Print the results ---
    print("--- Optimal Design Found ---")
    print(f"Container Type: {best_design['type']}")
    if best_design['type'] == 'Box':
        print(f"Dimensions (L x W x H): {best_design['dims'][0]}cm x {best_design['dims'][1]}cm x {best_design['dims'][2]}cm")
        print(f"Ball Arrangement: {' x '.join(map(str, best_design['layout']))}")
    else:
        print(f"Dimensions (Radius x Height): {best_design['dims'][0]}cm x {best_design['dims'][1]}cm")
        print(f"Ball Arrangement: {best_design['layout']}")

    print("\n--- Cost Calculation ---")
    print(f"Number of Energy Balls: {best_design['balls']}")
    print(f"Surface Area of Container: {best_design['area']:.2f} cm^2")
    print("\nTotal Cost = (Number of Balls * Cost per Ball) + (Surface Area * Cost per cm^2)")
    print(f"Total Cost = ({best_design['balls']} * {COST_PER_BALL_USD}) + ({best_design['area']:.2f} * {COST_PER_CM2_USD})")
    print(f"Total Cost = {best_design['cost_balls']} + {best_design['cost_container']:.2f}")
    print(f"Total Cost = ${best_design['cost']:.2f}")
    
    print(f"\nThe lowest total cost for the design (C) is: {final_cost}")
    print(f"<<<{final_cost}>>>")


solve_container_problem()
