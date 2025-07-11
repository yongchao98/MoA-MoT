import math

def solve_container_problem():
    """
    Solves the energy ball container design problem to find the minimum total cost.
    """
    # Step 1: Define constants
    ENERGY_PER_BALL_MJ = 25
    TOTAL_ENERGY_MJ = 1000
    COST_PER_BALL_USD = 1000
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = BALL_RADIUS_CM * 2
    MATERIAL_COST_PER_CM2 = 200
    PRECISION_CM = 0.5

    # Step 2: Calculate the number of balls and their fixed cost
    num_balls = math.ceil(TOTAL_ENERGY_MJ / ENERGY_PER_BALL_MJ)
    cost_balls = num_balls * COST_PER_BALL_USD

    print("--- Problem Analysis ---")
    print(f"Energy per ball: {ENERGY_PER_BALL_MJ} MJ")
    print(f"Required total energy: {TOTAL_ENERGY_MJ} MJ")
    print(f"Minimum number of balls required: {num_balls}")
    print(f"Fixed cost of energy balls: ${cost_balls:,.2f}")
    print("-" * 26)

    # --- Box Optimization ---
    min_box_area = float('inf')
    best_box_config = None
    # We search for arrangements of balls (nx, ny, nz)
    # The loop range is set slightly larger than num_balls to explore possibilities
    # where using more balls might create a more cube-like (efficient) shape.
    limit = int(num_balls * 1.5) 
    for nx in range(1, limit + 1):
        for ny in range(nx, limit + 1):
            if nx * ny > limit:
                break
            # Calculate the number of layers needed
            nz = math.ceil(num_balls / (nx * ny))
            
            l = nx * BALL_DIAMETER_CM
            w = ny * BALL_DIAMETER_CM
            h = nz * BALL_DIAMETER_CM
            
            area = 2 * (l*w + l*h + w*h)
            
            if area < min_box_area:
                min_box_area = area
                best_box_config = {
                    'type': 'Box',
                    'ball_config': (nx, ny, nz),
                    'total_balls': nx*ny*nz,
                    'dimensions': (l, w, h),
                    'area': area
                }

    # --- Cylinder Optimization ---
    min_cyl_area = float('inf')
    best_cyl_config = None
    # Search for arrangements (nx x ny base grid, nh layers high)
    for nx in range(1, limit + 1):
        for ny in range(nx, limit + 1):
            if nx * ny > limit:
                break
            # Calculate number of layers
            nh = math.ceil(num_balls / (nx * ny))

            # Ideal dimensions based on packing
            h_ideal = nh * BALL_DIAMETER_CM
            # Diameter must contain the diagonal of the nx*ny grid of balls
            d_ideal = BALL_DIAMETER_CM * math.sqrt(nx**2 + ny**2)

            # Apply manufacturing precision by rounding up to the nearest 0.5 cm
            h_final = math.ceil(h_ideal / PRECISION_CM) * PRECISION_CM
            d_final = math.ceil(d_ideal / PRECISION_CM) * PRECISION_CM
            r_final = d_final / 2.0
            
            area = (2 * math.pi * r_final**2) + (2 * math.pi * r_final * h_final)

            if area < min_cyl_area:
                min_cyl_area = area
                best_cyl_config = {
                    'type': 'Cylinder',
                    'ball_config': (nx, ny, nh),
                    'total_balls': nx*ny*nh,
                    'dimensions': (d_final, h_final),
                    'area': area
                }

    # --- Compare and Conclude ---
    cost_box_container = best_box_config['area'] * MATERIAL_COST_PER_CM2
    cost_cyl_container = best_cyl_config['area'] * MATERIAL_COST_PER_CM2
    
    total_cost_box = cost_balls + cost_box_container
    total_cost_cyl = cost_balls + cost_cyl_container
    
    print("\n--- Design Comparison ---")
    print(f"Best Box Design Cost: ${total_cost_box:,.2f}")
    print(f"Best Cylinder Design Cost: ${total_cost_cyl:,.2f}")
    print("-" * 26)

    if total_cost_box < total_cost_cyl:
        best_design = best_box_config
        container_cost = cost_box_container
        total_cost = total_cost_box
    else:
        best_design = best_cyl_config
        container_cost = cost_cyl_container
        total_cost = total_cost_cyl
        
    print("\n--- Optimal Design Chosen ---")
    print(f"The best container is a: {best_design['type']}")
    if best_design['type'] == 'Box':
        l, w, h = best_design['dimensions']
        nx, ny, nz = best_design['ball_config']
        print(f"Configuration: {nx}x{ny}x{nz} arrangement of balls")
        print(f"Dimensions (L x W x H): {l:.1f} cm x {w:.1f} cm x {h:.1f} cm")
    else: # Cylinder
        d, h = best_design['dimensions']
        nx, ny, nh = best_design['ball_config']
        print(f"Configuration: {nh} layers of a {nx}x{ny} ball grid")
        print(f"Dimensions (Diameter x Height): {d:.1f} cm x {h:.1f} cm")

    print(f"Container Surface Area: {best_design['area']:.2f} cm^2")
    
    print("\n--- Final Cost Calculation ---")
    print("Total Cost = Cost of Balls + Cost of Container")
    print(f"Total Cost = ({num_balls} balls * ${COST_PER_BALL_USD}/ball) + ({best_design['area']:.2f} cm^2 * ${MATERIAL_COST_PER_CM2}/cm^2)")
    print(f"Total Cost = ${cost_balls:,.2f} + ${container_cost:,.2f}")
    print(f"Total Cost = ${total_cost:,.2f}")
    
    return int(total_cost)

# Run the solver and print the final answer in the required format
final_cost = solve_container_problem()
print(f"<<<{final_cost}>>>")
