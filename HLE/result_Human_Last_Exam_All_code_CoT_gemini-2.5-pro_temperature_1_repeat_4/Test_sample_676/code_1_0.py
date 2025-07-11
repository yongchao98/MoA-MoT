import math

def solve_pioneer_probe_design():
    """
    Solves the energy ball container design problem to find the minimum total cost.
    """
    # --- Problem Constants ---
    ENERGY_PER_BALL = 25  # MJ
    COST_PER_BALL = 1000  # USD
    RADIUS_BALL = 2.0  # cm
    DIAMETER_BALL = 4.0 # cm
    COST_PER_CM2 = 200  # USD
    PRECISION = 0.5  # cm
    TARGET_ENERGY = 1000 # MJ

    min_total_cost = float('inf')
    best_config = {}

    # --- Step 1: Calculate minimum number of balls ---
    min_balls = math.ceil(TARGET_ENERGY / ENERGY_PER_BALL)

    # --- Step 6: Iterate through a range of ball counts ---
    # Start with the minimum required and check a few more.
    # Increasing balls adds significant cost, so a large increase is unlikely to be optimal.
    for n_balls in range(min_balls, min_balls + 25):
        cost_balls = n_balls * COST_PER_BALL

        # Optimization: if ball cost alone is higher than the current best total cost, stop.
        if cost_balls >= min_total_cost:
            break

        # --- Step 4: Optimize Box Container for n_balls ---
        min_box_area = float('inf')
        best_box_dims_balls = None
        
        # Search for the best grid (nx, ny, nz) to hold n_balls
        # To minimize surface area, dimensions should be as close as possible (cubic)
        for nx in range(1, n_balls + 1):
            # Optimization: if nx*1*1 is already too big, no need to check further
            if nx * 1 * 1 > n_balls and nx > n_balls**(1/3):
                break
            # To avoid duplicate permutations (e.g., 2x3x4 vs 3x2x4), ensure nx <= ny
            for ny in range(nx, n_balls + 1):
                # Optimization
                if nx * ny * 1 > n_balls and nx*ny > n_balls**(2/3):
                    break
                
                nz = math.ceil(n_balls / (nx * ny))
                
                area = 2 * ((nx * DIAMETER_BALL) * (ny * DIAMETER_BALL) + 
                              (nx * DIAMETER_BALL) * (nz * DIAMETER_BALL) + 
                              (ny * DIAMETER_BALL) * (nz * DIAMETER_BALL))
                
                if area < min_box_area:
                    min_box_area = area
                    best_box_dims_balls = (nx, ny, nz)
        
        total_cost_box = cost_balls + min_box_area * COST_PER_CM2
        if total_cost_box < min_total_cost:
            min_total_cost = total_cost_box
            L, W, H = (d * DIAMETER_BALL for d in best_box_dims_balls)
            best_config = {
                "type": "Box", "n_balls": n_balls, "cost_balls": cost_balls,
                "dims_cm": (L, W, H), "surface_area": min_box_area, "total_cost": total_cost_box,
                "final_equation": f"{cost_balls} + {min_box_area:.2f} * {COST_PER_CM2} = {total_cost_box:.2f}"
            }

        # --- Step 5: Optimize Cylinder Container for n_balls ---
        min_cyl_area = float('inf')
        best_cyl_dims_balls = None
        
        # Search for the best base grid (nx, ny) and number of layers (nz)
        for nx in range(1, n_balls + 1):
            if nx * 1 * 1 > n_balls and nx > n_balls**(1/3):
                break
            for ny in range(nx, n_balls + 1):
                if nx * ny * 1 > n_balls and nx*ny > n_balls**(2/3):
                    break

                nz = math.ceil(n_balls / (nx * ny))
                
                # The cylinder's base must enclose the rectangular grid of balls
                base_L = nx * DIAMETER_BALL
                base_W = ny * DIAMETER_BALL
                radius_unrounded = 0.5 * math.sqrt(base_L**2 + base_W**2)
                
                # Apply manufacturing precision
                radius = math.ceil(radius_unrounded / PRECISION) * PRECISION
                height = nz * DIAMETER_BALL
                
                area = 2 * math.pi * radius**2 + 2 * math.pi * radius * height
                
                if area < min_cyl_area:
                    min_cyl_area = area
                    best_cyl_dims_balls = (radius, height)

        total_cost_cyl = cost_balls + min_cyl_area * COST_PER_CM2
        if total_cost_cyl < min_total_cost:
            min_total_cost = total_cost_cyl
            best_config = {
                "type": "Cylinder", "n_balls": n_balls, "cost_balls": cost_balls,
                "dims_cm": best_cyl_dims_balls, "surface_area": min_cyl_area, "total_cost": total_cost_cyl,
                "final_equation": f"{cost_balls} + {min_cyl_area:.2f} * {COST_PER_CM2} = {total_cost_cyl:.2f}"
            }

    # --- Step 7: Output the final result ---
    if best_config:
        print(f"Optimal Design Found:")
        print(f"Container Type: {best_config['type']}")
        print(f"Number of Energy Balls: {best_config['n_balls']}")
        if best_config['type'] == 'Box':
            print(f"Container Dimensions (L x W x H): {best_config['dims_cm'][0]:.1f} cm x {best_config['dims_cm'][1]:.1f} cm x {best_config['dims_cm'][2]:.1f} cm")
        else: # Cylinder
            print(f"Container Dimensions (Radius x Height): {best_config['dims_cm'][0]:.1f} cm x {best_config['dims_cm'][1]:.1f} cm")
        print(f"Container Surface Area: {best_config['surface_area']:.2f} cm^2")
        print("\nCost Calculation:")
        print(f"Cost of Energy Balls = {best_config['n_balls']} * {COST_PER_BALL} = {best_config['cost_balls']}")
        material_cost = best_config['surface_area'] * COST_PER_CM2
        print(f"Cost of Container Material = {best_config['surface_area']:.2f} * {COST_PER_CM2} = {material_cost:.2f}")
        print(f"Total Cost = {best_config['cost_balls']} + {material_cost:.2f} = {best_config['total_cost']:.2f}")
        # The final answer C as requested
        print(f"\n<<<{int(round(best_config['total_cost']))}>>>")
    else:
        print("No solution found.")
        print("<<<0>>>")

solve_pioneer_probe_design()