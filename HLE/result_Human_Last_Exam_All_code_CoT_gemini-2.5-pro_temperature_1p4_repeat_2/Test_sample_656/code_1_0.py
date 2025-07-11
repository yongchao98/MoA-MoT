import math

def solve():
    """
    Solves the energy ball packing problem by searching for the optimal
    container dimensions for both box and cylinder shapes.
    """

    # --- Part 1: Box Optimization ---
    # We assume a simple cubic packing for large balls (radius 2cm, diameter 4cm).
    # To maximize the number of balls, dimensions should ideally be multiples of 4cm.
    # We scale everything by 2, so all dimensions are integers representing units of 0.5cm.
    # Large ball diameter = 4 cm = 8 units.
    # Surface area constraint: 2*(L*W + L*H + W*H) <= 1050 -> L_u*W_u + L_u*H_u + W_u*H_u <= 2100
    
    max_b_box = 0
    best_box_dims_u = None

    print("--- Searching for optimal box ---")
    # Iterate through dimensions in 0.5cm units (scaled by 2)
    # Assume L_u <= W_u <= H_u to avoid permutations
    # 3 * L_u^2 <= 2100 => L_u <= sqrt(700) approx 26.4. Scaled: 52.8
    for l_u in range(1, 53):
        # Rough bound for W_u: 2*W_u^2 < 2100 => W_u < 32. Scaled: 64
        for w_u in range(l_u, 65):
            if l_u * w_u > 2100:
                continue
            
            # Bound H_u based on the surface area constraint
            # h_u * (l_u + w_u) <= 2100 - l_u * w_u
            if (l_u + w_u) == 0:
                continue
            h_u_max = (2100 - l_u * w_u) / (l_u + w_u)
            if h_u_max < w_u:
                continue

            for h_u in range(w_u, int(h_u_max) + 1):
                # Number of large balls (diameter 8 units)
                b = (l_u // 8) * (w_u // 8) * (h_u // 8)

                if b > max_b_box:
                    max_b_box = b
                    best_box_dims_u = (l_u, w_u, h_u)
                    sa_term = l_u*w_u + l_u*h_u + w_u*h_u
                    print(f"New best box: dims_u={best_box_dims_u}, b={b}, sa={sa_term/2.0} cm^2")

    best_box_energy = max_b_box * 20

    # --- Part 2: Cylinder Optimization ---
    # We pack in layers (height 4cm). Each layer is a 2D circle packing problem.
    # Data for densest known packings of N circles of radius r in a circle of radius R.
    # Format: N -> required R/r ratio.
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cinc/cinc.html
    circle_packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.154, 4: 2.414, 5: 2.701, 6: 2.999, 7: 3.0, 
        8: 3.304, 9: 3.503, 10: 3.613, 11: 3.739, 12: 3.812, 13: 3.923, 
        14: 3.999, 15: 4.018, 16: 4.143, 17: 4.215, 18: 4.288, 19: 4.328, 
        20: 4.387, 21: 4.453, 22: 4.512, 23: 4.565, 24: 4.615, 25: 4.662,
        26: 4.708, 27: 4.742, 28: 4.819, 29: 4.863, 30: 4.918, 31: 4.954, 
        32: 5.0, 33: 5.039, 34: 5.088, 35: 5.121, 36: 5.161, 37: 5.215, 
        38: 5.25, 39: 5.291, 40: 5.323, 41: 5.353, 42: 5.393, 43: 5.437, 
        44: 5.479, 45: 5.499, 46: 5.549, 47: 5.587, 48: 5.626, 49: 5.659, 
        50: 5.698, 51: 5.742, 52: 5.772, 53: 5.823, 54: 5.861
    }
    
    # Invert for easier lookup: ratio -> max_n
    n_for_ratio = sorted(circle_packing_ratios.items(), key=lambda item: item[1])
    def get_n_for_ratio(ratio):
        max_n = 0
        for n, r_ratio in n_for_ratio:
            if r_ratio <= ratio:
                max_n = n
            else:
                break
        return max_n

    max_b_cyl = 0
    best_cyl_dims_u = None

    print("\n--- Searching for optimal cylinder ---")
    # Iterate through cylinder dimensions (scaled units)
    # SA constraint: R_u * (R_u + H_u) <= 2100/pi ~= 668.45
    for r_u in range(8, 51): # r_u >= 8 to fit one ball (diameter 8)
        # Calculate max h_u for this r_u
        h_u_max = (2100.0 / math.pi) / r_u - r_u
        if h_u_max < 8: # h_u must be at least 8 for one layer
            continue
        
        for h_u in range(8, int(h_u_max) + 1):
            num_layers = int((h_u / 2.0) / 4.0) # Height in cm / ball diameter
            if num_layers == 0:
                continue

            container_radius_cm = r_u / 2.0
            ball_radius_cm = 2.0
            ratio = container_radius_cm / ball_radius_cm
            
            n_per_layer = get_n_for_ratio(ratio)
            b = n_per_layer * num_layers

            if b > max_b_cyl:
                max_b_cyl = b
                best_cyl_dims_u = (r_u, h_u)
                sa = math.pi * (r_u / 2.0) * ( (r_u / 2.0) + (h_u / 2.0) )
                print(f"New best cylinder: dims_u={best_cyl_dims_u}, b={b}, sa={sa:.1f} cm^2")

    best_cyl_energy = max_b_cyl * 20

    # --- Part 3: Compare and finalize ---
    print("\n--- Comparison ---")
    print(f"Best box can hold {max_b_box} large balls for {best_box_energy} MJ.")
    print(f"Best cylinder can hold {max_b_cyl} large balls for {best_cyl_energy} MJ.")

    if best_cyl_energy > best_box_energy:
        r, h = best_cyl_dims_u[0] / 2.0, best_cyl_dims_u[1] / 2.0
        a = 0
        b = max_b_cyl
        final_energy = b * 20
        print("\nFinal Choice: Cylinder is better.")
        print(f"The container is a cylinder with r={r} cm, h={h} cm.")
        print(f"It can hold a={a} small balls and b={b} large balls.")
        print(f"Total potential energy packed: {a}*1 + {b}*20 = {final_energy} MJ")
        final_answer = f"[cylinder r={r}, h={h}]{a};{b}"

    else:
        l, w, h = best_box_dims_u[0]/2.0, best_box_dims_u[1]/2.0, best_box_dims_u[2]/2.0
        a = 0
        b = max_b_box
        final_energy = b * 20
        print("\nFinal Choice: Box is better.")
        print(f"The container is a box with dimensions {l}x{w}x{h} cm.")
        print(f"It can hold a={a} small balls and b={b} large balls.")
        print(f"Total potential energy packed: {a}*1 + {b}*20 = {final_energy} MJ")
        final_answer = f"[box {l}x{w}x{h}]{a};{b}"

    print("\n<<<" + final_answer + ">>>")


solve()