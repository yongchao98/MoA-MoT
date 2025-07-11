import math

def solve_container_problem():
    """
    Solves the energy ball container optimization problem.
    """
    # --- Constants ---
    TOTAL_ENERGY_REQ = 1000  # MJ
    ENERGY_PER_BALL = 25     # MJ
    COST_PER_BALL = 1000     # USD
    RADIUS_BALL = 2          # cm
    DIAMETER_BALL = 4        # cm
    COST_PER_AREA = 200      # USD per cm^2
    PRECISION = 0.5          # cm

    # --- Step 1: Calculate Ball Requirements ---
    num_balls = math.ceil(TOTAL_ENERGY_REQ / ENERGY_PER_BALL)
    cost_balls = num_balls * COST_PER_BALL

    print("Step 1: Calculate the number of energy balls and their cost.")
    print(f"Total energy required: {TOTAL_ENERGY_REQ} MJ")
    print(f"Energy per ball: {ENERGY_PER_BALL} MJ")
    print(f"Number of balls required = ceil({TOTAL_ENERGY_REQ} / {ENERGY_PER_BALL}) = {num_balls}")
    print(f"Cost of energy balls = {num_balls} balls * ${COST_PER_BALL}/ball = ${cost_balls:,.2f}")
    print("\n" + "="*40 + "\n")

    # --- Step 2: Find the optimal 'Box' container design ---
    # By analyzing factorizations of 40, the most cube-like arrangement is 5x4x2.
    # This minimizes surface area for a given volume.
    n_l, n_w, n_h = 5, 4, 2
    L = n_l * DIAMETER_BALL
    W = n_w * DIAMETER_BALL
    H = n_h * DIAMETER_BALL
    area_box = 2 * (L*W + W*H + H*L)
    cost_box_container = area_box * COST_PER_AREA
    total_cost_box = cost_balls + cost_box_container

    print("Step 2: Find the optimal 'Box' container design.")
    print(f"The optimal grid arrangement to pack {num_balls} balls is {n_l} x {n_w} x {n_h}.")
    print(f"This requires a box of dimensions (L x W x H): {L:.1f} cm x {W:.1f} cm x {H:.1f} cm.")
    print(f"Surface Area = 2 * ({L:.1f}*{W:.1f} + {W:.1f}*{H:.1f} + {H:.1f}*{L:.1f}) = {area_box:.2f} cm^2")
    print(f"Box container cost = {area_box:.2f} cm^2 * ${COST_PER_AREA}/cm^2 = ${cost_box_container:,.2f}")
    print(f"Total cost with box design = ${cost_balls:,.2f} + ${cost_box_container:,.2f} = ${total_cost_box:,.2f}")
    print("\n" + "="*40 + "\n")

    # --- Step 3: Find the optimal 'Cylinder' container design ---
    # We test arrangements with k balls per layer.
    # We use known minimal radius ratios (R_container / r_ball) for packing k circles.
    # Ratios for k=1 to 12: {k: R/r}
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0,
        7: 3.0, 8: 3.305, 9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029
    }
    
    min_area_cyl = float('inf')
    best_cyl_config = {}

    for k in range(1, num_balls + 1):
        if k in packing_ratios:
            # Calculate required radius and round up to precision
            min_r_container = packing_ratios[k] * RADIUS_BALL
            R_cyl = math.ceil(min_r_container / PRECISION) * PRECISION
            
            # Calculate required height
            num_layers = math.ceil(num_balls / k)
            H_cyl = num_layers * DIAMETER_BALL
            
            area_cyl = 2 * math.pi * R_cyl**2 + 2 * math.pi * R_cyl * H_cyl
            
            if area_cyl < min_area_cyl:
                min_area_cyl = area_cyl
                best_cyl_config = {
                    'k': k, 'R': R_cyl, 'H': H_cyl, 'area': area_cyl, 'layers': num_layers
                }

    area_cyl = best_cyl_config['area']
    cost_cyl_container = area_cyl * COST_PER_AREA
    total_cost_cyl = cost_balls + cost_cyl_container

    print("Step 3: Find the optimal 'Cylinder' container design.")
    print(f"Optimal arrangement found by testing layers of k=1 to {num_balls} balls.")
    print(f"The best design has {best_cyl_config['layers']} layers with up to {best_cyl_config['k']} balls per layer.")
    print(f"This requires a cylinder of Radius {best_cyl_config['R']:.1f} cm and Height {best_cyl_config['H']:.1f} cm.")
    print(f"Surface Area = 2*pi*R^2 + 2*pi*R*H = 2*pi*{best_cyl_config['R']:.1f}^2 + 2*pi*{best_cyl_config['R']:.1f}*{best_cyl_config['H']:.1f} = {area_cyl:.2f} cm^2")
    print(f"Cylinder container cost = {area_cyl:.2f} cm^2 * ${COST_PER_AREA}/cm^2 = ${cost_cyl_container:,.2f}")
    print(f"Total cost with cylinder design = ${cost_balls:,.2f} + ${cost_cyl_container:,.2f} = ${total_cost_cyl:,.2f}")
    print("\n" + "="*40 + "\n")

    # --- Step 4: Compare and Conclude ---
    print("Step 4: Compare costs and determine the final design.")
    if total_cost_cyl < total_cost_box:
        best_design_name = "Cylinder"
        min_total_cost = total_cost_cyl
        final_area = area_cyl
        final_R = best_cyl_config['R']
        final_H = best_cyl_config['H']
        print(f"Comparing the two designs, the Cylinder is cheaper by ${total_cost_box - total_cost_cyl:,.2f}.")
        print(f"The minimum total cost C is ${min_total_cost:,.2f}")
        print("\nThe final cost C is calculated using the Cylinder's parameters:")
        print(f"C = (Number of Balls * Cost per Ball) + (Container Surface Area * Cost per cm^2)")
        print(f"C = ({num_balls} * {COST_PER_BALL}) + ({cost_per_area} * (2 * pi * {final_R:.1f}^2 + 2 * pi * {final_R:.1f} * {final_H:.1f}))")
    else:
        best_design_name = "Box"
        min_total_cost = total_cost_box
        final_area = area_box
        print(f"Comparing the two designs, the Box is cheaper by ${total_cost_cyl - total_cost_box:,.2f}.")
        print(f"The minimum total cost C is ${min_total_cost:,.2f}")
        print("\nThe final cost C is calculated using the Box's parameters:")
        print(f"C = (Number of Balls * Cost per Ball) + (Container Surface Area * Cost per cm^2)")
        print(f"C = ({num_balls} * {COST_PER_BALL}) + ({cost_per_area} * (2 * ({L:.1f}*{W:.1f} + {W:.1f}*{H:.1f} + {H:.1f}*{L:.1f})))")

    # Final Answer
    final_C = min_total_cost
    print(f"\n<<<C = {final_C:.2f}>>>")

solve_container_problem()