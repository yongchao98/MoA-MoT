import math

def solve_tomb_engraving():
    """
    Calculates the optimal number of squares (N) and circles (M) to maximize
    the total number of engraved characters (K) on a 140x110cm metal sheet.
    """

    # --- Step 1: Define Constants ---
    mat_w, mat_h = 140, 110
    sq_w = 10
    sq_chars = 4
    circ_r = 20
    circ_d = 40

    # --- Step 2: Calculate Characters per Circle ---
    # To encode 441 unique characters with 8-state Bagua symbols, we need:
    # ceil(log8(441)) = 3 symbols per character.
    symbols_per_char = math.ceil(math.log(441, 8))
    # A circle holds 999 symbols.
    circ_chars = math.floor(999 / symbols_per_char)

    print("Step 1: Determine the number of characters per artifact.")
    print(f" - Each square (10x10cm) holds {sq_chars} characters.")
    print(f" - Each circle (20cm radius) holds {circ_chars} characters.\n")

    # --- Step 3 & 4: Iterate through configurations to solve the packing problem ---
    print("Step 2: Analyze packing configurations to maximize K = (N * 4) + (M * 333).")

    # We pre-calculate bounding boxes for M circles using compact packing arrangements.
    # Format: (M, bbox_width, bbox_height, description)
    configs = [
        (0, 0, 0, "0 circles"),
        (1, circ_d, circ_d, "a single circle"),
        (2, 2 * circ_d, circ_d, "2 circles in a line"),
        (3, 3 * circ_d, circ_d, "3 circles in a line"),
        (4, 2 * circ_d, 2 * circ_d, "a 2x2 grid of circles"),
        (5, 3 * circ_d, 2 * circ_d, "a cluster of 5 fitting in a 120x80 box"),
        (6, 6 * circ_r, (2 + math.sqrt(3)) * circ_r, "a hexagonal pack of 6 circles"),
        (7, 6 * circ_r, 2 * (1 + math.sqrt(3)) * circ_r, "a hexagonal pack of 7 circles")
    ]
    
    best_N, best_M, best_K = -1, -1, -1

    for M, bbox_w, bbox_h, desc in configs:
        # Check if the bounding box fits on the material. Since 140>110 and all our
        # bboxes have one dimension <= 110 or not, we only need one orientation check.
        if bbox_w > mat_w or bbox_h > mat_h:
            print(f" - For M={M}, the bounding box ({bbox_w:.2f}x{bbox_h:.2f}cm) does not fit.")
            continue

        # Place the circle cluster's bounding box in a corner of the material.
        # This leaves two main rectangular areas for squares.
        # Rectangle 1: Adjacent to the bbox's width
        rem1_w, rem1_h = mat_w - bbox_w, mat_h
        n1 = math.floor(rem1_w / sq_w) * math.floor(rem1_h / sq_w)

        # Rectangle 2: Above the bbox
        rem2_w, rem2_h = bbox_w, mat_h - bbox_h
        n2 = math.floor(rem2_w / sq_w) * math.floor(rem2_h / sq_w)
        
        N = n1 + n2
        K = (N * sq_chars) + (M * circ_chars)

        print(f" - For M={M} ({desc}), the bbox is {bbox_w:.2f}x{bbox_h:.2f}cm. "
              f"This leaves space for N={N} squares. Total K = {K}.")

        if K > best_K:
            best_K, best_N, best_M = K, N, M
            
    # --- Step 5: Output the Final Result ---
    print("\nStep 3: Identify the best configuration.")
    print(f"The configuration that maximizes the total number of characters is M={best_M} circles and N={best_N} squares.")
    
    final_equation = f"The maximal number of characters is ({best_N} * {sq_chars}) + ({best_M} * {circ_chars}) = {best_N * sq_chars} + {best_M * circ_chars} = {best_K}."
    print(final_equation)

    final_answer = f"{best_N}:{best_M}:{best_K}"
    print(f"\nThe final answer in the format N:M:K is: {final_answer}")

solve_tomb_engraving()
<<<22:7:2419>>>