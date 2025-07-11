import math

def solve_pioneer_packing():
    """
    This function solves the container design problem for the Pioneer probe.
    The strategy is to find the optimal box dimensions to maximize the packed energy,
    based on the principle of prioritizing the more energy-dense large balls.
    """

    # --- Constants based on the problem description ---
    BALL_A_DIAMETER = 4.0
    BALL_A_ENERGY = 20
    MAX_SURFACE_AREA = 1050.0

    # --- Step 1: Optimize for large balls in a box container ---
    # To maximize the number of large balls (4cm diameter), we assume the box
    # dimensions L, W, H are multiples of 4. So, L=4*k_l, W=4*k_w, H=4*k_h.
    # The surface area constraint 2*(LW+LH+WH) <= 1050 becomes:
    # 32 * (k_l*k_w + k_l*k_h + k_w*k_h) <= 1050
    # or k_l*k_w + k_l*k_h + k_w*k_h <= 32.8125
    MAX_SA_COEFF = MAX_SURFACE_AREA / 32.0

    best_k_product = 0
    best_k = (0, 0, 0)
    
    # We search for the integer triplet (k_l, k_w, k_h) that maximizes the
    # product k_l*k_w*k_h (proportional to the number of balls) while
    # satisfying the surface area constraint.
    limit = 15 # A reasonable search limit for k values
    for k_l in range(1, limit):
        for k_w in range(k_l, limit):
            for k_h in range(k_w, limit):
                sa_coeff = k_l*k_w + k_l*k_h + k_w*k_h
                if sa_coeff <= MAX_SA_COEFF:
                    k_product = k_l * k_w * k_h
                    if k_product > best_k_product:
                        best_k_product = k_product
                        best_k = (k_l, k_w, k_h)

    # --- Step 2: Calculate container dimensions and ball count ---
    k_l, k_w, k_h = best_k
    
    # Optimal dimensions for the box
    L = int(BALL_A_DIAMETER * k_l)
    W = int(BALL_A_DIAMETER * k_w)
    H = int(BALL_A_DIAMETER * k_h)
    
    # Number of large balls (b, 2-cm radius)
    num_large_balls = k_l * k_w * k_h
    
    # Since the dimensions are perfect multiples of the large ball diameter,
    # there is no leftover space to fit any small balls.
    num_small_balls = 0
    
    # --- Step 3: Print the results and the final formatted answer ---
    print("--- Optimal Container Design ---")
    print(f"Container Type: box")
    print(f"Dimensions (L x W x H): {L}x{W}x{H} cm")
    print(f"Surface Area: {2 * (L*W + L*H + W*H)} cm^2 (Limit: {MAX_SURFACE_AREA} cm^2)")
    print(f"Number of 1-cm radius balls (a): {num_small_balls}")
    print(f"Number of 2-cm radius balls (b): {num_large_balls}")
    print(f"Maximum Potential Energy: {num_large_balls * BALL_A_ENERGY} MJ")
    
    container_desc = f"box {L}x{W}x{H}"
    a = num_small_balls
    b = num_large_balls
    final_answer = f"[{container_desc}]{a};{b}"
    
    print("\n--- Formatted Answer ---")
    print(final_answer)

# Execute the solver
solve_pioneer_packing()