import math

def solve_container_problem():
    """
    Calculates the minimum cost to build a container for energy balls,
    considering both box and cylinder shapes.
    """

    # --- Problem Constants ---
    TOTAL_ENERGY_MJ = 1000
    BALL_ENERGY_MJ = 25
    BALL_COST_USD = 1000
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    MATERIAL_COST_PER_CM2 = 200
    PRECISION_CM = 0.5
    PI = math.pi

    # --- Minimum Balls Calculation ---
    min_balls_required = math.ceil(TOTAL_ENERGY_MJ / BALL_ENERGY_MJ)

    # --- Box Optimization ---
    best_box = {'cost': float('inf')}
    
    # For a box, the surface area is minimized for a fixed number of balls.
    # Increasing the number of balls will always increase the minimum possible surface area.
    # Therefore, we only need to check for the minimum required number of balls.
    num_balls_box = min_balls_required
    min_surface_area_box = float('inf')
    best_factors = None

    for i in range(1, int(num_balls_box**0.5) + 1):
        if num_balls_box % i == 0:
            for j in range(i, int((num_balls_box/i)**0.5) + 1):
                if (num_balls_box / i) % j == 0:
                    k = num_balls_box // (i * j)
                    factors = (i, j, k)
                    l = factors[0] * BALL_DIAMETER_CM
                    w = factors[1] * BALL_DIAMETER_CM
                    h = factors[2] * BALL_DIAMETER_CM
                    area = 2 * (l*w + w*h + h*l)
                    if area < min_surface_area_box:
                        min_surface_area_box = area
                        best_factors = factors

    cost_balls_box = num_balls_box * BALL_COST_USD
    cost_material_box = min_surface_area_box * MATERIAL_COST_PER_CM2
    total_cost_box = cost_balls_box + cost_material_box
    best_box = {
        'cost': total_cost_box,
        'type': 'Box',
        'balls': num_balls_box,
        'dims': (best_factors[0] * BALL_DIAMETER_CM, best_factors[1] * BALL_DIAMETER_CM, best_factors[2] * BALL_DIAMETER_CM),
        'area': min_surface_area_box,
        'cost_material': cost_material_box,
        'cost_balls': cost_balls_box,
        'ball_arrangement': best_factors
    }

    # --- Cylinder Optimization ---
    best_cylinder = {'cost': float('inf')}
    
    # Radii of smallest circle enclosing k unit circles (Rc). Source: Specht, E. (2021) "The best known packings of equal circles in a circle"
    # We need R_layer = BALL_RADIUS_CM * Rc
    packing_radii_Rc = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 4.864, 20: 5.122
    }

    for k, Rc in packing_radii_Rc.items():
        # Calculate cylinder dimensions based on packing k balls per layer
        required_radius = BALL_RADIUS_CM * Rc
        
        # Enforce 0.5 cm precision for radius
        cyl_radius = math.ceil(required_radius / PRECISION_CM) * PRECISION_CM
        
        num_layers = math.ceil(min_balls_required / k)
        cyl_height = num_layers * BALL_DIAMETER_CM
        
        # Calculate total balls and costs
        num_balls_cyl = k * num_layers
        surface_area_cyl = (2 * PI * cyl_radius**2) + (2 * PI * cyl_radius * cyl_height)
        
        cost_balls_cyl = num_balls_cyl * BALL_COST_USD
        cost_material_cyl = surface_area_cyl * MATERIAL_COST_PER_CM2
        total_cost_cyl = cost_balls_cyl + cost_material_cyl
        
        if total_cost_cyl < best_cylinder['cost']:
            best_cylinder = {
                'cost': total_cost_cyl,
                'type': 'Cylinder',
                'balls': num_balls_cyl,
                'dims': (cyl_radius, cyl_height),
                'area': surface_area_cyl,
                'cost_material': cost_material_cyl,
                'cost_balls': cost_balls_cyl,
                'k_per_layer': k,
                'num_layers': num_layers
            }

    # --- Compare and Print Final Result ---
    if best_box['cost'] < best_cylinder['cost']:
        winner = best_box
    else:
        winner = best_cylinder

    print(f"The optimal design is a {winner['type']}.")
    print("-" * 30)

    if winner['type'] == 'Box':
        print(f"Ball Arrangement (l x w x h): {winner['ball_arrangement'][0]} x {winner['ball_arrangement'][1]} x {winner['ball_arrangement'][2]}")
        print(f"Container Dimensions (L x W x H): {winner['dims'][0]} cm x {winner['dims'][1]} cm x {winner['dims'][2]} cm")
        print(f"Container Surface Area: {winner['area']:.2f} cm^2")
        print(f"Number of Energy Balls: {winner['balls']}")
        print("\nCost Calculation:")
        print(f"Material Cost = {winner['area']:.2f} cm^2 * ${MATERIAL_COST_PER_CM2}/cm^2 = ${winner['cost_material']:.2f}")
        print(f"Energy Ball Cost = {winner['balls']} balls * ${BALL_COST_USD}/ball = ${winner['cost_balls']:.2f}")
        print(f"Total Cost = ${winner['cost_material']:.2f} + ${winner['cost_balls']:.2f} = ${winner['cost']:.2f}")
    
    else: # Cylinder
        print(f"Balls per Layer: {winner['k_per_layer']}")
        print(f"Number of Layers: {winner['num_layers']}")
        print(f"Container Dimensions (Radius x Height): {winner['dims'][0]:.1f} cm x {winner['dims'][1]:.1f} cm")
        print(f"Container Surface Area: {winner['area']:.4f} cm^2")
        print(f"Number of Energy Balls: {winner['balls']}")
        print("\nCost Calculation:")
        print(f"Material Cost = {winner['area']:.4f} cm^2 * ${MATERIAL_COST_PER_CM2}/cm^2 = ${winner['cost_material']:.4f}")
        print(f"Energy Ball Cost = {winner['balls']} balls * ${BALL_COST_USD}/ball = ${winner['cost_balls']:.2f}")
        print(f"Total Cost = ${winner['cost_material']:.4f} + ${winner['cost_balls']:.2f} = ${winner['cost']:.4f}")
        
    final_cost = winner['cost']
    print(f"\n<<<C = {final_cost}>>>")


solve_container_problem()