import math

def solve_container_problem():
    """
    Solves the container design optimization problem by calculating the minimum cost
    for both box and cylinder designs and choosing the best one.
    """

    # --- Constants ---
    BALL_DIAMETER_CM = 4.0
    BALL_RADIUS_CM = 2.0
    BALL_ENERGY_MJ = 25
    BALL_COST_USD = 1000
    MATERIAL_COST_PER_CM2_USD = 200
    REQUIRED_ENERGY_MJ = 1000
    PRECISION_CM = 0.5

    # --- Step 1: Calculate the number of balls and their fixed cost ---
    num_balls = int(math.ceil(REQUIRED_ENERGY_MJ / BALL_ENERGY_MJ))
    cost_balls = num_balls * BALL_COST_USD

    # --- Step 2: Find the best Box design ---
    def find_best_box_design(n_balls):
        """Finds the box with the minimum surface area for n_balls."""
        best_area = float('inf')
        
        # Find integer factorizations of n_balls: i x j x k
        for i in range(1, n_balls + 1):
            if n_balls % i == 0:
                for j in range(i, n_balls + 1):
                    if (n_balls % (i * j)) == 0:
                        k = n_balls // (i * j)
                        if k >= j:
                            # We have factors (i, j, k) for ball counts
                            L = i * BALL_DIAMETER_CM
                            W = j * BALL_DIAMETER_CM
                            H = k * BALL_DIAMETER_CM
                            
                            area = 2 * (L * W + W * H + H * L)
                            if area < best_area:
                                best_area = area
        return best_area

    # --- Step 3: Find the best Cylinder design ---
    def find_best_cylinder_design(n_balls):
        """Finds the cylinder with the minimum surface area for n_balls."""

        # Known optimal packing ratios (container radius / circle radius) for k circles
        # Source: packomania.com and other mathematical sources for circle packing
        ratios = {
            1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
            8: 3.305, 9: 3.613, 10: 3.813, 11: 3.924, 12: 4.030, 13: 4.236,
            14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 4.864,
            20: 5.122, 21: 5.323, 22: 5.348, 23: 5.517, 24: 5.606, 25: 5.680, 
            26: 5.759, 27: 5.864, 28: 5.992, 29: 6.094, 30: 6.104, 31: 6.223, 
            32: 6.281, 33: 6.388, 34: 6.402, 35: 6.474, 36: 6.478, 37: 6.478, 
            38: 6.690, 39: 6.746, 40: 6.799
        }

        best_area = float('inf')

        # Iterate through possible numbers of layers
        for n_layers in range(1, n_balls + 1):
            k_per_layer = math.ceil(n_balls / n_layers)
            
            if k_per_layer > n_balls: continue

            # Calculate container radius and apply precision
            ratio = ratios.get(k_per_layer, math.sqrt(k_per_layer / 0.9069)) # Use approx for >40
            min_radius = ratio * BALL_RADIUS_CM
            container_radius = math.ceil(min_radius / PRECISION_CM) * PRECISION_CM
            
            # Height is based on simple stacking, which meets precision requirement
            container_height = n_layers * BALL_DIAMETER_CM

            area = 2 * math.pi * container_radius * (container_radius + container_height)
            
            if area < best_area:
                best_area = area

        return best_area

    # --- Step 4: Compare designs and print the result ---
    box_area = find_best_box_design(num_balls)
    cyl_area = find_best_cylinder_design(num_balls)

    box_material_cost = box_area * MATERIAL_COST_PER_CM2_USD
    cyl_material_cost = cyl_area * MATERIAL_COST_PER_CM2_USD
    
    box_total_cost = box_material_cost + cost_balls
    cyl_total_cost = cyl_material_cost + cost_balls
    
    print("Design Analysis:")
    print("-" * 20)
    print(f"Number of energy balls needed: {num_balls}")
    print(f"Cost of energy balls: {num_balls} balls * ${BALL_COST_USD}/ball = ${cost_balls:,.2f}")
    print("-" * 20)
    print("Box Container Option:")
    print(f"Minimum surface area: {box_area:,.2f} cm^2")
    print(f"Total cost: ${box_area:,.2f} * ${MATERIAL_COST_PER_CM2_USD} + ${cost_balls:,.2f} = ${box_total_cost:,.2f}")
    print("-" * 20)
    print("Cylinder Container Option:")
    print(f"Minimum surface area: {cyl_area:,.2f} cm^2")
    print(f"Total cost: ${cyl_area:,.2f} * ${MATERIAL_COST_PER_CM2_USD} + ${cost_balls:,.2f} = ${cyl_total_cost:,.2f}")
    print("-" * 20)
    
    # Final Choice
    if box_total_cost < cyl_total_cost:
        design_choice = "Box"
        final_cost = box_total_cost
        final_area = box_area
        material_cost = box_material_cost
    else:
        design_choice = "Cylinder"
        final_cost = cyl_total_cost
        final_area = cyl_area
        material_cost = cyl_material_cost

    print(f"\nThe optimal design is the {design_choice} container.")
    print("The final cost C is calculated as C = (Surface Area * Material Cost) + (Number of Balls * Ball Cost)")
    print(f"C = ({final_area:.2f} cm^2 * ${MATERIAL_COST_PER_CM2_USD}/cm^2) + ({num_balls} balls * ${BALL_COST_USD}/ball)")
    print(f"C = ${material_cost:,.2f} + ${cost_balls:,.2f}")
    print(f"C = ${final_cost:,.2f}")
    
    # Final answer in the required format
    print(f"\n<<<{final_cost:.2f}>>>")

# Execute the solver
solve_container_problem()
