import math

def solve_container_problem():
    """
    Solves the container optimization problem by calculating the minimum cost for
    both a box and a cylinder container, then selecting the cheaper option.
    """

    # --- Problem Parameters ---
    ENERGY_PER_BALL_MJ = 25
    TOTAL_ENERGY_TARGET_MJ = 1000
    COST_PER_BALL_USD = 1000
    MATERIAL_COST_PER_CM2_USD = 200
    BALL_RADIUS_CM = 2
    BALL_DIAMETER_CM = 4
    PRECISION_CM = 0.5

    # --- Step 1: Calculate Ball Requirements ---
    num_balls = math.ceil(TOTAL_ENERGY_TARGET_MJ / ENERGY_PER_BALL_MJ)
    cost_balls = num_balls * COST_PER_BALL_USD

    # --- Function to calculate Box Cost ---
    def calculate_box_cost():
        min_surface_area = float('inf')
        best_arrangement = (0, 0, 0)

        # Iterate through possible arrangements (nx, ny, nz) of balls
        # We assume packing exactly num_balls is optimal to avoid extra ball costs
        limit = int(num_balls**(1/3.0)) + 2 # Heuristic limit for search
        for nx in range(1, num_balls + 1):
            if nx > limit * 2: continue # Optimization
            for ny in range(nx, num_balls + 1):
                if nx * ny > num_balls:
                    break
                nz = math.ceil(num_balls / (nx * ny))
                if nz < ny:
                    continue

                # Box dimensions determined by ball diameter
                L = nx * BALL_DIAMETER_CM
                W = ny * BALL_DIAMETER_CM
                H = nz * BALL_DIAMETER_CM
                
                # Check precision (always a multiple of 4cm, so it's a multiple of 0.5cm)
                surface_area = 2 * (L * W + L * H + W * H)
                
                if surface_area < min_surface_area:
                    min_surface_area = surface_area
                    best_arrangement = (L, W, H)

        cost_container = min_surface_area * MATERIAL_COST_PER_CM2_USD
        total_cost = cost_container + cost_balls
        return total_cost, best_arrangement, min_surface_area

    # --- Function to calculate Cylinder Cost ---
    def calculate_cylinder_cost():
        min_surface_area = float('inf')
        best_config = {}

        # Radii of minimal enclosing circle for n circles of radius 1
        # Data from packomania.com for n=1 to 40
        optimal_radii_r1 = [
            0.0, 1.0, 2.0, 2.155, 2.414, 2.701, 3.0, 3.0, 3.305, 3.507, 3.813, 
            3.924, 3.963, 4.030, 4.289, 4.328, 4.521, 4.615, 4.693, 4.864, 4.864, 
            5.122, 5.176, 5.303, 5.340, 5.517, 5.632, 5.687, 5.864, 5.986, 6.083, 
            6.149, 6.208, 6.208, 6.436, 6.467, 6.467, 6.467, 6.727, 6.747, 6.941
        ]

        # Iterate through possible number of layers, k
        for k in range(1, num_balls + 1):
            n_layer = math.ceil(num_balls / k)
            if n_layer > len(optimal_radii_r1) - 1:
                continue

            # Calculate required dimensions
            radius_ratio = optimal_radii_r1[n_layer]
            required_diameter = BALL_RADIUS_CM * 2 * radius_ratio
            height = k * BALL_DIAMETER_CM

            # Adjust for manufacturing precision
            final_diameter = math.ceil(required_diameter / PRECISION_CM) * PRECISION_CM
            final_radius = final_diameter / 2
            
            # Height is a multiple of 4, so it meets precision
            final_height = height

            surface_area = 2 * math.pi * final_radius**2 + 2 * math.pi * final_radius * final_height
            
            if surface_area < min_surface_area:
                min_surface_area = surface_area
                best_config = {
                    "radius": final_radius, 
                    "height": final_height,
                    "layers": k,
                    "balls_per_layer": n_layer
                }
        
        cost_container = min_surface_area * MATERIAL_COST_PER_CM2_USD
        total_cost = cost_container + cost_balls
        return total_cost, best_config, min_surface_area

    # --- Main Logic: Compare designs and print result ---
    box_cost, box_dims, box_area = calculate_box_cost()
    cyl_cost, cyl_config, cyl_area = calculate_cylinder_cost()

    if box_cost < cyl_cost:
        print("The optimal design is a Box.")
        print(f"Number of balls (N) = {num_balls}")
        print(f"Ball arrangement (nx*ny*nz) requires a container of size:")
        print(f"Box Dimensions (L, W, H) = {box_dims[0]:.1f} cm x {box_dims[1]:.1f} cm x {box_dims[2]:.1f} cm")
        print(f"Container Surface Area = {box_area:.2f} cm^2")
        print(f"Cost of Balls = {num_balls} * ${COST_PER_BALL_USD:.2f} = ${cost_balls:.2f}")
        container_cost = box_area * MATERIAL_COST_PER_CM2_USD
        print(f"Cost of Container Material = {box_area:.2f} cm^2 * ${MATERIAL_COST_PER_CM2_USD:.2f}/cm^2 = ${container_cost:.2f}")
        print(f"Total Cost = Cost_Balls + Cost_Container = ${cost_balls:.2f} + ${container_cost:.2f} = ${box_cost:.2f}")
        final_answer = box_cost
    else:
        print("The optimal design is a Cylinder.")
        print(f"Number of balls (N) = {num_balls}")
        print(f"Configuration: {cyl_config['layers']} layers, with up to {cyl_config['balls_per_layer']} balls per layer.")
        print(f"Cylinder Height (H) = {cyl_config['height']:.1f} cm")
        print(f"Cylinder Radius (R) = {cyl_config['radius']:.2f} cm")
        print(f"Container Surface Area = {cyl_area:.2f} cm^2")
        print(f"Cost of Balls = {num_balls} * ${COST_PER_BALL_USD:.2f} = ${cost_balls:.2f}")
        container_cost = cyl_area * MATERIAL_COST_PER_CM2_USD
        print(f"Cost of Container Material = {cyl_area:.2f} cm^2 * ${MATERIAL_COST_PER_CM2_USD:.2f}/cm^2 = ${container_cost:.2f}")
        print(f"Total Cost = Cost_Balls + Cost_Container = ${cost_balls:.2f} + ${container_cost:.2f} = ${cyl_cost:.2f}")
        final_answer = cyl_cost

    print(f"\n<<<{final_answer:.2f}>>>")

if __name__ == '__main__':
    solve_container_problem()
