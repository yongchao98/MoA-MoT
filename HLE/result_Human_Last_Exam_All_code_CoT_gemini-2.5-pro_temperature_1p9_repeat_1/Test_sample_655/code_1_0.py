import math

def solve_packing():
    """
    Finds the optimal container design and ball packing to maximize energy.

    This function searches for the optimal box dimensions that stay within
    the surface area constraint of 1050 cm^2. It assumes a cubic packing
    model for the spherical energy balls (4x4x4 for large, 2x2x2 for small)
    to make the problem computationally feasible. A greedy algorithm is used,
    prioritizing the more energy-dense large balls.
    """
    best_energy = 0
    best_config = {}

    # Iterate through possible dimensions (L, W, H) in increments of 0.5 cm.
    # To reduce search space, we enforce L >= W >= H and set a practical upper limit.
    # We iterate using integers and divide by 2.0 to handle 0.5 cm precision.
    # A limit of 50 corresponds to a max dimension of 25.0 cm.
    for l_int in range(1, 51):
        l = l_int / 2.0
        for w_int in range(1, l_int + 1):
            w = w_int / 2.0
            for h_int in range(1, w_int + 1):
                h = h_int / 2.0

                # Constraint 1: Check surface area
                surface_area = 2 * (l * w + l * h + w * h)
                if surface_area > 1050:
                    continue

                # --- Start Greedy Packing Calculation ---

                # Step 1: Pack as many large balls (4x4x4 cubes) as possible.
                n2_l = math.floor(l / 4)
                n2_w = math.floor(w / 4)
                n2_h = math.floor(h / 4)
                n2 = n2_l * n2_w * n2_h

                # Dimensions of the block occupied by large balls
                l_used_by_n2 = n2_l * 4
                w_used_by_n2 = n2_w * 4
                h_used_by_n2 = n2_h * 4

                # Step 2: Pack small balls (2x2x2 cubes) in the remaining space.
                # The remaining space is decomposed into 3 non-overlapping slabs.
                n1 = 0

                # Slab 1: Along the length dimension
                rem_l = l - l_used_by_n2
                if rem_l >= 2:
                    n1 += math.floor(rem_l / 2) * math.floor(w / 2) * math.floor(h / 2)
                
                # Slab 2: Along the width dimension
                rem_w = w - w_used_by_n2
                if rem_w >= 2:
                    n1 += math.floor(l_used_by_n2 / 2) * math.floor(rem_w / 2) * math.floor(h / 2)

                # Slab 3: Along the height dimension
                rem_h = h - h_used_by_n2
                if rem_h >= 2:
                    n1 += math.floor(l_used_by_n2 / 2) * math.floor(w_used_by_n2 / 2) * math.floor(rem_h / 2)

                # --- End Greedy Packing Calculation ---

                # Calculate total energy and update the best configuration if needed
                current_energy = 10 * n2 + n1
                if current_energy > best_energy:
                    best_energy = current_energy
                    best_config = {
                        "shape_desc": f"box {l}x{w}x{h}",
                        "n1": n1,
                        "n2": n2,
                    }

    # Print the final result in the specified format
    if best_config:
        # [C]a;b with a is the number of 1-cm balls, b is the number of 2-cm balls
        print(f"{best_config['shape_desc']};{best_config['n1']};{best_config['n2']}")
    else:
        print("No solution found within the given constraints.")

solve_packing()
<<<box 14.0x12.5x12.5;36;27>>>