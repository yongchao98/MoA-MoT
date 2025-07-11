import math

def solve_container_problem():
    """
    Solves the container optimization problem by analyzing box and cylinder designs.
    """
    # 1. Define constants and initial calculations
    ENERGY_NEEDED = 1000  # MJ
    ENERGY_PER_BALL = 25  # MJ
    COST_PER_BALL = 1000  # usd
    BALL_RADIUS = 2.0     # cm
    BALL_DIAMETER = 4.0   # cm
    MATERIAL_COST_PER_CM2 = 200 # usd
    PRECISION = 0.5       # cm

    # Calculate the number of balls and their fixed cost
    NUM_BALLS = math.ceil(ENERGY_NEEDED / ENERGY_PER_BALL)
    BALL_COST = NUM_BALLS * COST_PER_BALL

    # --- 2. Optimize Box Container ---
    min_box_sa = float('inf')
    best_box_dims_cm = (0, 0, 0)
    best_box_counts = (0, 0, 0)

    # Search for the best integer grid (l, w, h) to arrange the balls
    # We search for l, w, h such that l*w*h >= NUM_BALLS
    # To minimize surface area, l, w, h should be close to the cube root of NUM_BALLS
    limit = int(NUM_BALLS**(1/3.0)) + 3 
    for l in range(1, NUM_BALLS + 1):
        for w in range(l, NUM_BALLS + 1):
            if l * w > NUM_BALLS * 2: # Optimization: avoid very flat shapes
                break
            h = math.ceil(NUM_BALLS / (l * w))
            
            L, W, H = l * BALL_DIAMETER, w * BALL_DIAMETER, h * BALL_DIAMETER
            sa = 2 * (L * W + W * H + H * L)
            
            if sa < min_box_sa:
                min_box_sa = sa
                best_box_dims_cm = (L, W, H)
                best_box_counts = (l,w,h)

    box_container_cost = min_box_sa * MATERIAL_COST_PER_CM2
    total_box_cost = BALL_COST + box_container_cost
    
    # --- 3. Optimize Cylinder Container ---
    
    # R/r ratios for optimal packing of N circles in a circle (N -> R/r).
    # This data is from known solutions to the circle packing problem.
    # R = BALL_RADIUS * ratio
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0, 8: 3.307,
        9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236, 14: 4.328, 15: 4.521,
        16: 4.641, 17: 4.747, 18: 4.864, 19: 4.868, 20: 5.122, 21: 5.209, 22: 5.312,
        23: 5.344, 24: 5.46, 25: 5.603, 26: 5.617, 27: 5.688, 28: 5.801, 29: 5.86,
        30: 5.922, 31: 6.035, 32: 6.095, 33: 6.148, 34: 6.257, 35: 6.307, 36: 6.407,
        37: 6.471, 38: 6.516, 39: 6.561, 40: 6.666
    }
    
    def ceil_to_precision(value, precision):
        return math.ceil(value / precision) * precision

    min_cyl_sa = float('inf')
    best_cyl_dims = (0, 0)

    # Iterate through possible number of layers to stack the balls
    for layers in range(1, NUM_BALLS + 1):
        balls_per_layer = math.ceil(NUM_BALLS / layers)
        if balls_per_layer > max(packing_ratios.keys()):
            continue

        H_container = layers * BALL_DIAMETER
        
        ratio = packing_ratios[balls_per_layer]
        R_ideal = BALL_RADIUS * ratio
        R_container = ceil_to_precision(R_ideal, PRECISION)

        sa = 2 * math.pi * R_container**2 + 2 * math.pi * R_container * H_container
        
        if sa < min_cyl_sa:
            min_cyl_sa = sa
            best_cyl_dims = (R_container, H_container)

    cyl_container_cost = min_cyl_sa * MATERIAL_COST_PER_CM2
    total_cyl_cost = BALL_COST + cyl_container_cost

    # --- 4. Compare and Print Final Result ---
    if total_cyl_cost < total_box_cost:
        final_design = "Cylinder"
        final_cost = total_cyl_cost
        radius, height = best_cyl_dims
        surface_area = min_cyl_sa
        container_cost = cyl_container_cost
        
        print(f"The best design is a {final_design} container.")
        print(f"Number of balls: {NUM_BALLS}")
        print(f"Cylinder Radius: {radius:.1f} cm")
        print(f"Cylinder Height: {height:.1f} cm")
        print(f"Cylinder Surface Area: {surface_area:.2f} cm^2")
        print("\n--- Cost Calculation ---")
        print("Total Cost = (Number of balls * Cost per ball) + (Surface Area * Material Cost per cm^2)")
        print(f"Total Cost = ({NUM_BALLS} * {COST_PER_BALL}) + ({surface_area:.2f} * {MATERIAL_COST_PER_CM2})")
        print(f"Total Cost = {BALL_COST} + {container_cost:.2f}")
        
    else:
        final_design = "Box"
        final_cost = total_box_cost
        length, width, height = sorted(best_box_dims_cm, reverse=True)
        surface_area = min_box_sa
        container_cost = box_container_cost

        print(f"The best design is a {final_design} container.")
        print(f"Number of balls: {NUM_BALLS}")
        print(f"Box Dimensions: {length:.1f} cm x {width:.1f} cm x {height:.1f} cm")
        print(f"Box Surface Area: {surface_area:.2f} cm^2")
        print("\n--- Cost Calculation ---")
        print("Total Cost = (Number of balls * Cost per ball) + (Surface Area * Material Cost per cm^2)")
        print(f"Total Cost = ({NUM_BALLS} * {COST_PER_BALL}) + ({surface_area:.2f} * {MATERIAL_COST_PER_CM2})")
        print(f"Total Cost = {BALL_COST} + {container_cost:.2f}")

    final_cost_int = round(final_cost)
    print(f"Final Total Cost = {final_cost_int}")
    print(f"<<<{final_cost_int}>>>")


solve_container_problem()
