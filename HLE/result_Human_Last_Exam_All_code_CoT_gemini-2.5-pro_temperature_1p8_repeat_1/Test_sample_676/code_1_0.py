import math

def solve_pioneer_container_problem():
    """
    Solves the container design optimization problem for the Pioneer probe.
    """
    # 1. Problem Parameters
    BALL_RADIUS_CM = 2.0
    BALL_DIAMETER_CM = 2 * BALL_RADIUS_CM
    BALL_ENERGY_MJ = 25
    BALL_COST_USD = 1000.0
    TOTAL_ENERGY_REQ_MJ = 1000.0
    MATERIAL_COST_PER_CM2 = 200.0
    PRECISION_CM = 0.5

    # 2. Calculate fixed ball cost
    num_balls = math.ceil(TOTAL_ENERGY_REQ_MJ / BALL_ENERGY_MJ)
    total_ball_cost = num_balls * BALL_COST_USD

    # --- Box Optimization ---
    best_box_area = float('inf')
    best_box_dims = None

    # Possible grid arrangements (nx, ny, nz) where nx*ny*nz = 40
    arrangements = [
        (40, 1, 1), (20, 2, 1), (10, 4, 1),
        (10, 2, 2), (8, 5, 1), (5, 4, 2)
    ]

    for nx, ny, nz in arrangements:
        L = nx * BALL_DIAMETER_CM
        W = ny * BALL_DIAMETER_CM
        H = nz * BALL_DIAMETER_CM
        area = 2 * (L * W + L * H + W * H)
        if area < best_box_area:
            best_box_area = area
            best_box_dims = (L, W, H)

    box_material_cost = best_box_area * MATERIAL_COST_PER_CM2

    # --- Cylinder Optimization ---
    best_cyl_area = float('inf')
    best_cyl_dims = None
    best_cyl_params = None
    
    # Radii of smallest enclosing circle for k unit circles (r=1)
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cinc/cinc.html
    # Values R(k) for packing k circles of radius 1
    packing_radii = {
        1: 1.0, 2: 2.0, 3: 2.1547, 4: 2.4142, 5: 2.7013, 6: 3.0, 7: 3.0, 
        8: 3.3048, 9: 3.6131, 10: 3.8130, 11: 3.9238, 12: 4.0296, 13: 4.2361, 
        14: 4.3284, 15: 4.5213, 16: 4.6154, 17: 4.7923, 18: 4.8637, 19: 5.0,
        20: 5.1223, 21: 5.2393, 22: 5.3183, 23: 5.5036, 24: 5.5948, 25: 5.6881, 
        26: 5.7629, 27: 5.8637, 28: 6.0, 29: 6.0772, 30: 6.1437, 31: 6.2361, 
        32: 6.3572, 33: 6.4566, 34: 6.5168, 35: 6.6433, 36: 6.7475, 37: 6.7475, 
        38: 6.8797, 39: 7.0, 40: 7.0583
    }

    for k in range(1, int(num_balls) + 1):
        # Calculate radius required for k balls in a layer
        r_unit_packing = packing_radii[k]
        required_radius = r_unit_packing * BALL_RADIUS_CM
        
        # Apply precision: round up to the nearest multiple of PRECISION_CM
        container_radius = math.ceil(required_radius / PRECISION_CM) * PRECISION_CM
        
        # Calculate height required
        num_layers = math.ceil(num_balls / k)
        container_height = num_layers * BALL_DIAMETER_CM
        
        area = (2 * math.pi * container_radius * container_height) + (2 * math.pi * container_radius**2)
        
        if area < best_cyl_area:
            best_cyl_area = area
            best_cyl_dims = (container_radius, container_height)
            best_cyl_params = k

    cyl_material_cost = best_cyl_area * MATERIAL_COST_PER_CM2

    # 4. Compare and find the absolute best design
    if box_material_cost < cyl_material_cost:
        print("Optimal Design: Box Container")
        L, W, H = best_box_dims
        final_material_cost = box_material_cost
        final_area = best_box_area
        print(f"Dimensions (L, W, H): ({L:.1f} cm, {W:.1f} cm, {H:.1f} cm)")
        print(f"Container Surface Area = 2 * ({L:.1f}*{W:.1f} + {L:.1f}*{H:.1f} + {W:.1f}*{H:.1f}) = {final_area:.4f} cm^2")
    else:
        print("Optimal Design: Cylinder Container")
        R, H = best_cyl_dims
        final_material_cost = cyl_material_cost
        final_area = best_cyl_area
        print(f"Dimensions (Radius, Height): ({R:.1f} cm, {H:.1f} cm)")
        print(f"Configuration: {math.ceil(num_balls / best_cyl_params)} layers of up to {best_cyl_params} balls each.")
        print(f"Container Surface Area = 2 * pi * {R:.1f} * {H:.1f} + 2 * pi * {R:.1f}^2 = {final_area:.4f} cm^2")

    total_cost = total_ball_cost + final_material_cost
    
    print("\n--- Cost Calculation ---")
    print(f"Number of Energy Balls: {int(num_balls)}")
    print(f"Material Cost = {final_area:.4f} cm^2 * ${MATERIAL_COST_PER_CM2:.2f}/cm^2 = ${final_material_cost:.2f}")
    print(f"Ball Cost = {int(num_balls)} balls * ${BALL_COST_USD:.2f}/ball = ${total_ball_cost:.2f}")
    print("\n--- Final Equation ---")
    print(f"Total Cost = Ball Cost + Material Cost")
    print(f"Total Cost = ${total_ball_cost:.2f} + ${final_material_cost:.2f} = ${total_cost:.2f}")
    
    # Final answer format
    print(f"\n<<<{total_cost:.2f}>>>")

# Execute the solver
solve_pioneer_container_problem()
