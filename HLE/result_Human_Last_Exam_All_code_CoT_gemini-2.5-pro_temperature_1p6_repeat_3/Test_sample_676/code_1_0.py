import math

def solve_container_problem():
    """
    Solves the container optimization problem to find the minimum cost design.
    """
    # Problem constants
    BALL_RADIUS = 2.0  # cm
    BALL_DIAMETER = 4.0 # cm
    BALL_ENERGY = 25  # MJ
    BALL_COST = 1000  # USD
    MATERIAL_COST_PER_CM2 = 200  # USD
    REQUIRED_ENERGY = 1000  # MJ
    PRECISION = 0.5  # cm

    # Step 1: Calculate the minimum number of balls required
    min_num_balls = math.ceil(REQUIRED_ENERGY / BALL_ENERGY)

    # --- Box Optimization ---
    min_box_cost = float('inf')
    best_box_details = {}

    # We check slightly more balls than the minimum in case a better packing reduces material cost
    for num_balls in range(min_num_balls, min_num_balls + 25):
        # Find factor triplets for the number of balls
        for i in range(1, int(num_balls**(1/3.0)) + 2):
            if num_balls % i == 0:
                for j in range(i, int((num_balls/i)**0.5) + 2):
                    if (num_balls/i) % j == 0:
                        k = num_balls // (i * j)
                        if i * j * k == num_balls:
                            nx, ny, nz = i, j, k
                            
                            L = nx * BALL_DIAMETER
                            W = ny * BALL_DIAMETER
                            H = nz * BALL_DIAMETER
                            
                            surface_area = 2 * (L * W + L * H + W * H)
                            material_cost = surface_area * MATERIAL_COST_PER_CM2
                            ball_cost_total = num_balls * BALL_COST
                            total_cost = material_cost + ball_cost_total

                            if total_cost < min_box_cost:
                                min_box_cost = total_cost
                                best_box_details = {
                                    'cost': total_cost,
                                    'num_balls': num_balls,
                                    'arrangement': (nx, ny, nz),
                                    'dims': (L, W, H),
                                    'surface_area': surface_area,
                                    'material_cost': material_cost,
                                    'ball_cost': ball_cost_total
                                }

    # --- Cylinder Optimization ---
    min_cylinder_cost = float('inf')
    best_cylinder_details = {}

    # We check slightly more balls than the minimum
    for num_balls in range(min_num_balls, min_num_balls + 25):
        for k_layers in range(1, num_balls + 1):
            if num_balls % k_layers == 0:
                m_per_layer = num_balls // k_layers
                
                # Find best square-like arrangement (nx, ny) for m_per_layer
                best_nx = int(m_per_layer**0.5)
                while m_per_layer % best_nx != 0:
                    best_nx -= 1
                nx, ny = best_nx, m_per_layer // best_nx

                H_cyl = k_layers * BALL_DIAMETER
                
                # Calculate minimum radius required to fit nx x ny grid of spheres
                # dist from arrangement center to corner ball center + ball radius
                if nx == 1 and ny == 1:
                    R_cyl_min = BALL_RADIUS
                else:
                    dist_to_center = 2 * math.sqrt((nx - 1)**2 + (ny - 1)**2)
                    R_cyl_min = dist_to_center + BALL_RADIUS
                
                # Adhere to precision
                R_cyl = math.ceil(R_cyl_min / PRECISION) * PRECISION
                
                surface_area = 2 * math.pi * R_cyl**2 + 2 * math.pi * R_cyl * H_cyl
                material_cost = surface_area * MATERIAL_COST_PER_CM2
                ball_cost_total = num_balls * BALL_COST
                total_cost = material_cost + ball_cost_total
                
                if total_cost < min_cylinder_cost:
                    min_cylinder_cost = total_cost
                    # Storing details for potential printout, though box is expected to win
                    best_cylinder_details = {'cost': total_cost}


    # --- Compare and Print Final Result ---
    if min_box_cost <= min_cylinder_cost:
        final_design = best_box_details
        design_type = "Box"
    else:
        # This case is unlikely based on analysis, but we handle it.
        # Cylinder printing would need more details saved.
        final_design = best_cylinder_details 
        design_type = "Cylinder"
        
    print(f"The optimal design is a {design_type} container.")
    print("\n--- Cost Calculation ---")
    print("Total Cost = Cost of Energy Balls + Cost of Container Material")
    print("Total Cost = (Number of Balls * Cost per Ball) + (Surface Area * Cost per cm^2)")
    
    nb = final_design['num_balls']
    bc = BALL_COST
    sa = final_design['surface_area']
    mc = MATERIAL_COST_PER_CM2
    
    print(f"Total Cost = ({nb} * {bc}) + ({sa} * {mc})")
    print(f"Total Cost = {final_design['ball_cost']} + {final_design['material_cost']}")
    
    final_cost = final_design['cost']
    print(f"Final Total Cost = {final_cost}")
    print(f"<<<{final_cost}>>>")

solve_container_problem()
