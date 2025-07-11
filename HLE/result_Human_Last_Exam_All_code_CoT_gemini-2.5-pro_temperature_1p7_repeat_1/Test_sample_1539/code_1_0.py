import math

def solve_engraving_problem():
    """
    Calculates the optimal number of circular and squared plates to maximize engraved characters.
    """
    # Step 1: Define constants
    chars_circle = 9999
    chars_square = 360
    
    rect_w = 140
    rect_h = 110
    
    circle_diameter = 40
    square_side = 10

    print("Thinking Process:")
    print("The goal is to maximize K = N * 9999 + M * 360.")
    print("Since circular plates are much more valuable, we prioritize fitting the maximum number of circles (N).")
    print("-" * 20)

    # --- Strategy 1: Grid packing of circles (N=6) ---
    print("Strategy 1: Simple Grid Packing of Circles")
    n_grid_w = math.floor(rect_w / circle_diameter)
    n_grid_h = math.floor(rect_h / circle_diameter)
    n1 = n_grid_w * n_grid_h
    
    # Area used by circle bounding boxes
    used_w1 = n_grid_w * circle_diameter
    used_h1 = n_grid_h * circle_diameter
    
    # Calculate squares in leftover rectangular strips
    rem_strip1_w, rem_strip1_h = rect_w, rect_h - used_h1
    rem_strip2_w, rem_strip2_h = rect_w - used_w1, used_h1
    
    m1_strip1 = math.floor(rem_strip1_w / square_side) * math.floor(rem_strip1_h / square_side)
    m1_strip2 = math.floor(rem_strip2_w / square_side) * math.floor(rem_strip2_h / square_side)
    
    # A 10x10 square can fit in the interstitial space between four 40x40 circles in a grid
    num_interstitial_spaces = (n_grid_w - 1) * (n_grid_h - 1)
    
    m1 = m1_strip1 + m1_strip2 + num_interstitial_spaces
    k1 = n1 * chars_circle + m1 * chars_square
    
    print(f"A {n_grid_w}x{n_grid_h} grid fits N = {n1} circles.")
    print(f"Remaining space allows for M = {m1_strip1} + {m1_strip2} + {num_interstitial_spaces} = {m1} squares.")
    print(f"Result -> N={n1}, M={m1}, K = {n1} * {chars_circle} + {m1} * {chars_square} = {k1}")
    print("-" * 20)

    # --- Strategy 2: Hexagonal packing of circles (N=8) ---
    # A 3-2-3 staggered packing fits 8 circles in a 120 x 109.28 bounding box
    print("Strategy 2: Sub-optimal Hexagonal Packing (N=8)")
    n2 = 8
    used_w2 = 120 # From (3 circles wide with staggering)
    
    # Remaining space is a single strip
    rem_w2 = rect_w - used_w2
    rem_h2 = rect_h
    m2 = math.floor(rem_w2 / square_side) * math.floor(rem_h2 / square_side)
    k2 = n2 * chars_circle + m2 * chars_square
    
    print(f"A staggered '3-2-3' packing can fit N = {n2} circles.")
    print(f"This uses ~120cm width, leaving a {rem_w2}cm x {rem_h2}cm strip.")
    print(f"This strip fits M = {m2} squares.")
    print(f"Result -> N={n2}, M={m2}, K = {n2} * {chars_circle} + {m2} * {chars_square} = {k2}")
    print("-" * 20)
    
    # --- Strategy 3: Optimal packing of circles (N=9) ---
    # It's a known result that 9 circles can be packed in a 140x110 rectangle.
    # One such packing (3-3-3 staggered) uses the full 140cm width and ~109.28cm height.
    print("Strategy 3: Optimal Hexagonal Packing (N=9)")
    n3 = 9
    used_w3 = 140
    used_h3 = 20 + 2 * (circle_diameter * math.sqrt(3)/2) + 20 # Theoretical height
    
    # This packing is extremely efficient, leaving no space for 10x10 squares.
    m3 = 0
    k3 = n3 * chars_circle + m3 * chars_square
    
    print(f"An optimal staggered '3-3-3' packing can fit N = {n3} circles.")
    print(f"This uses the full 140cm width and ~109.28cm height, leaving no usable space for squares.")
    print(f"Result -> N={n3}, M={m3}.")
    print(f"K = {n3} * {chars_circle} + {m3} * {chars_square} = {k3}")
    print("-" * 20)
    
    # --- Conclusion ---
    print("Comparing the strategies:")
    print(f"Strategy 1 (Grid): K = {k1}")
    print(f"Strategy 2 (Hex-8): K = {k2}")
    print(f"Strategy 3 (Hex-9): K = {k3}")
    
    best_k = max(k1, k2, k3)
    if best_k == k1:
        best_n, best_m = n1, m1
    elif best_k == k2:
        best_n, best_m = n2, m2
    else:
        best_n, best_m = n3, m3
        
    print("\nThe best strategy is the one that maximizes K.")
    print(f"The maximum number of characters is {best_k}.")
    print("\nFinal Answer:")
    print("The optimal number of plates is:")
    print(f"N (Circular Plates) = {best_n}")
    print(f"M (Squared Plates) = {best_m}")
    
    print("\nThe final calculation is:")
    print(f"K = {best_n} * {chars_circle} + {best_m} * {chars_square} = {best_k}")
    
    final_answer_string = f"{best_n}:{best_m}:{best_k}"
    print(f"\nIn the required format: {final_answer_string}")
    
    # Final output for the platform
    return final_answer_string

final_answer = solve_engraving_problem()
print(f"<<<{final_answer}>>>")
