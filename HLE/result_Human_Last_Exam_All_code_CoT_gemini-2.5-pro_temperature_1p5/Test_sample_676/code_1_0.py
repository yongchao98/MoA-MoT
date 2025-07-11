import math

def solve_container_problem():
    """
    Solves the energy ball container optimization problem.
    """

    # --- Problem Constants ---
    TOTAL_ENERGY_REQ_MJ = 1000
    ENERGY_PER_BALL_MJ = 25
    COST_PER_BALL_USD = 1000
    COST_PER_AREA_USD_CM2 = 200
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    PRECISION_CM = 0.5

    # --- Step 1: Determine the number of balls ---
    num_balls = math.ceil(TOTAL_ENERGY_REQ_MJ / ENERGY_PER_BALL_MJ)

    def optimize_box(n_balls):
        """Finds the optimal box dimensions for n_balls."""
        min_area = float('inf')
        best_dims = (0, 0, 0)
        best_layout = (0, 0, 0)

        # Iterate through possible layouts (nx, ny, nz)
        # Search is pruned to be more efficient
        for nx in range(1, n_balls + 1):
            if nx * 1 * 1 > n_balls: break
            for ny in range(nx, n_balls + 1):
                if nx * ny * 1 > n_balls: break
                # Smallest nz that fits the balls
                nz = math.ceil(n_balls / (nx * ny))
                
                length = nx * BALL_DIAMETER_CM
                width = ny * BALL_DIAMETER_CM
                height = nz * BALL_DIAMETER_CM
                
                # Box dimensions are automatically multiples of 0.5cm since diameter is 4cm
                
                area = 2 * (length * width + length * height + width * height)
                
                if area < min_area:
                    min_area = area
                    best_dims = (length, width, height)
                    best_layout = (nx, ny, nz)
        
        return min_area, best_dims, best_layout

    def optimize_cylinder(n_balls):
        """Finds the optimal cylinder dimensions for n_balls."""
        min_area = float('inf')
        best_dims = (0, 0)
        best_k = 0

        # Packing ratios (R/r) for k circles in a circle. R=container radius, r=circle radius.
        # Source: http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html
        packing_ratios = {
            1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
            8: 3.304, 9: 3.503, 10: 3.813, 11: 3.924, 12: 4.030, 13: 4.236,
            14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 4.864,
            20: 5.122, 21: 5.210, 22: 5.347, 23: 5.514, 24: 5.611, 25: 5.753,
            26: 5.867, 27: 5.962, 28: 6.136, 29: 6.225, 30: 6.326, 31: 6.458,
            32: 6.545, 33: 6.669, 34: 6.748, 35: 6.840, 36: 6.993, 37: 7.001,
            38: 7.146, 39: 7.234, 40: 7.332
        }

        for k in range(1, n_balls + 1):
            num_layers = math.ceil(n_balls / k)
            
            # Calculate raw height and apply precision
            raw_height = num_layers * BALL_DIAMETER_CM
            # Height is automatically a multiple of 0.5cm
            height = raw_height
            
            # Calculate raw radius and apply precision
            ratio = packing_ratios.get(k, k) # Fallback for k > 40
            raw_radius = ratio * BALL_RADIUS_CM
            radius = math.ceil(raw_radius / PRECISION_CM) * PRECISION_CM
            
            area = (2 * math.pi * radius**2) + (2 * math.pi * radius * height)

            if area < min_area:
                min_area = area
                best_dims = (radius, height)
                best_k = k
        
        return min_area, best_dims, best_k

    # --- Step 2: Perform Optimizations ---
    box_area, box_dims, box_layout = optimize_box(num_balls)
    cyl_area, cyl_dims, cyl_k = optimize_cylinder(num_balls)

    # --- Step 3: Compare results and calculate total cost ---
    if box_area < cyl_area:
        design_type = "Box"
        min_surface_area = box_area
        final_dimensions = box_dims
    else:
        design_type = "Cylinder"
        min_surface_area = cyl_area
        final_dimensions = cyl_dims

    cost_balls = num_balls * COST_PER_BALL_USD
    cost_container = min_surface_area * COST_PER_AREA_USD_CM2
    total_cost = cost_balls + cost_container

    # --- Step 4: Print the results ---
    print("--- Design Optimization Report ---")
    print(f"Minimum number of energy balls required for >= {TOTAL_ENERGY_REQ_MJ} MJ: {num_balls}")
    print("\n--- Comparing Container Designs ---")
    print(f"Best Box Design (Layout {box_layout}):")
    print(f"  - Dimensions (LxWxH): {box_dims[0]:.1f} x {box_dims[1]:.1f} x {box_dims[2]:.1f} cm")
    print(f"  - Surface Area: {box_area:.2f} cm^2")
    print(f"Best Cylinder Design (Balls/layer k={cyl_k}):")
    print(f"  - Dimensions (RxH): {cyl_dims[0]:.1f} x {cyl_dims[1]:.1f} cm")
    print(f"  - Surface Area: {cyl_area:.2f} cm^2")
    
    print("\n--- Optimal Design Choice ---")
    print(f"The best design is a {design_type} container.")
    if design_type == "Cylinder":
      print(f"Radius = {final_dimensions[0]:.1f} cm, Height = {final_dimensions[1]:.1f} cm")
    else:
      print(f"Length = {final_dimensions[0]:.1f} cm, Width = {final_dimensions[1]:.1f} cm, Height = {final_dimensions[2]:.1f} cm")

    print(f"Minimum Surface Area = {min_surface_area:.4f} cm^2")
    
    print("\n--- Total Cost Calculation ---")
    print(f"Cost of Balls = {num_balls} balls * ${COST_PER_BALL_USD}/ball = ${cost_balls:.2f}")
    print(f"Cost of Container = {min_surface_area:.4f} cm^2 * ${COST_PER_AREA_USD_CM2}/cm^2 = ${cost_container:.2f}")
    
    print("\n--- Final Cost Equation ---")
    print(f"Total Cost = {num_balls} * {COST_PER_BALL_USD} + {min_surface_area} * {COST_PER_AREA_USD_CM2}")

    print(f"\nFinal Total Cost = ${total_cost:.2f}")
    
    # Required for submission format
    global C
    C = total_cost

# Execute the solver
solve_container_problem()
<<<266194.67>>>