import math

def solve_container_problem():
    """
    Solves the energy ball container design problem by calculating the minimum cost
    for both a box and a cylinder container, and choosing the cheaper option.
    """

    # --- Problem Constants ---
    ENERGY_PER_BALL_MJ = 25
    REQUIRED_ENERGY_MJ = 1000
    BALL_COST_USD = 1000
    MATERIAL_COST_USD_PER_CM2 = 200
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    PRECISION_CM = 0.5

    # Step 1: Calculate the minimum number of balls required
    num_balls = math.ceil(REQUIRED_ENERGY_MJ / ENERGY_PER_BALL_MJ)
    ball_total_cost = num_balls * BALL_COST_USD

    print(f"Minimum number of energy balls required: {num_balls}")
    print(f"Total cost of balls: ${ball_total_cost}\n")

    # --- Step 2: Box Container Optimization ---
    print("--- Analyzing Box Container ---")
    
    # We need to pack 40 balls. To minimize surface area, we need dimensions to be as close as possible.
    # We find integer triplets (nx,ny,nz) such that nx*ny*nz = 40.
    # The factors of 40 are 1, 2, 4, 5, 8, 10, 20, 40.
    # Triplets: (1,1,40), (1,2,20), (1,4,10), (1,5,8), (2,2,10), (2,4,5).
    # The most cube-like triplet is (2,4,5).
    nx, ny, nz = 2, 4, 5
    
    # Dimensions of the box
    box_l = nz * BALL_DIAMETER_CM
    box_w = ny * BALL_DIAMETER_CM
    box_h = nx * BALL_DIAMETER_CM
    
    # Surface Area of the box
    sa_box = 2 * (box_l * box_w + box_w * box_h + box_h * box_l)
    
    # Total cost for the box design
    material_cost_box = sa_box * MATERIAL_COST_USD_PER_CM2
    total_cost_box = ball_total_cost + material_cost_box
    
    print(f"Best arrangement: {nx}x{ny}x{nz} balls")
    print(f"Box dimensions (cm): L={box_l}, W={box_w}, H={box_h}")
    print(f"Box surface area: {sa_box:.2f} cm^2")
    print(f"Box material cost: ${material_cost_box:.2f}")
    print(f"Total cost for Box Design: ${total_cost_box:.2f}\n")

    # --- Step 3: Cylinder Container Optimization ---
    print("--- Analyzing Cylinder Container ---")

    min_cyl_sa = float('inf')
    best_cyl_config = {}

    # Ratios of container radius to inner circle radius for optimal 2D packing
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cinc/cinc.html
    k_ratios = {
        1: 1.0, 2: 2.0, 3: 1 + 2/math.sqrt(3), 4: 1 + math.sqrt(2),
        5: 2.7013, 6: 3.0, 7: 3.0, 8: 3.3047, 9: 3.6131, 10: 3.8130,
        11: 3.9237, 12: 4.029, 13: 4.236, 14: 4.328, 15: 4.521
    }

    # Iterate through possible numbers of balls per layer (k)
    for k in range(1, num_balls + 1):
        if k not in k_ratios:
            continue
            
        num_layers = math.ceil(num_balls / k)
        
        # Calculate radius and apply precision
        ideal_r = k_ratios[k] * BALL_RADIUS_CM
        cyl_r = math.ceil(ideal_r / PRECISION_CM) * PRECISION_CM
        
        # Calculate height based on stacking type and apply precision
        # Use denser (HCP-like) stacking for hexagonal arrangements (k=3,7)
        if k == 7 or k == 3:
            # Vertical distance between layers in dense packing
            layer_height = BALL_DIAMETER_CM * math.sqrt(2/3)
            ideal_h = (num_layers - 1) * layer_height + BALL_DIAMETER_CM
        else: # simple stacking
            ideal_h = num_layers * BALL_DIAMETER_CM

        cyl_h = math.ceil(ideal_h / PRECISION_CM) * PRECISION_CM
        
        sa_cyl = 2 * math.pi * cyl_r**2 + 2 * math.pi * cyl_r * cyl_h
        
        if sa_cyl < min_cyl_sa:
            min_cyl_sa = sa_cyl
            best_cyl_config = {
                'k': k, 'layers': num_layers, 'R': cyl_r, 'H': cyl_h
            }
    
    sa_cylinder = min_cyl_sa
    material_cost_cylinder = sa_cylinder * MATERIAL_COST_USD_PER_CM2
    total_cost_cylinder = ball_total_cost + material_cost_cylinder

    print(f"Best arrangement: {best_cyl_config['layers']} layers of {best_cyl_config['k']} balls/layer capacity")
    print(f"Cylinder dimensions (cm): R={best_cyl_config['R']}, H={best_cyl_config['H']}")
    print(f"Cylinder surface area: {sa_cylinder:.2f} cm^2")
    print(f"Cylinder material cost: ${material_cost_cylinder:.2f}")
    print(f"Total cost for Cylinder Design: ${total_cost_cylinder:.2f}\n")

    # --- Step 4: Final Comparison ---
    print("--- Final Comparison ---")
    if total_cost_cylinder < total_cost_box:
        final_cost = total_cost_cylinder
        cost_type = "Cylinder"
        final_sa = sa_cylinder
    else:
        final_cost = total_cost_box
        cost_type = "Box"
        final_sa = sa_box

    final_cost_int = round(final_cost)
    
    print(f"The best design is the {cost_type} container.")
    print("\nFinal Cost Calculation:")
    print(f"Total Cost = (Number of Balls * Cost per Ball) + (Surface Area * Material Cost per cm^2)")
    # The user asked to see each number in the final equation
    print(f"Total Cost = ({num_balls} * ${BALL_COST_USD}) + ({final_sa:.2f} cm^2 * ${MATERIAL_COST_USD_PER_CM2}) = ${final_cost:.2f}")
    print(f"The minimum total cost, C, is ${final_cost_int}.")
    
    # Final answer format
    print(f"\n<<<{final_cost_int}>>>")

solve_container_problem()