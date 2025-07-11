import math

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize the
    number of engraved characters, and prints the step-by-step solution.
    """

    # --- Problem Parameters ---
    mat_w = 140  # Material width in cm
    mat_h = 110  # Material height in cm

    sq_size = 10         # Square side length in cm
    chars_per_sq = 4     # Characters on one square

    circ_diameter = 40   # Circle diameter in cm (bounding box size)
    chars_per_circ = 999 # Characters on one circle

    # --- Analysis ---
    print("This is a packing problem to maximize the number of engraved characters, K.")
    print(f"The objective is to maximize K = ({chars_per_sq} * N) + ({chars_per_circ} * M).")
    print("Given the high value of characters per circle, we prioritize maximizing M (circles).")
    print("-" * 40)

    # --- Step 1: Calculate the maximum number of circles (M) ---
    print("Step 1: Calculate the maximum number of circles (M)")
    # We fit the circle's 40x40cm bounding box into the 140x110cm material.
    circles_along_w = math.floor(mat_w / circ_diameter)
    circles_along_h = math.floor(mat_h / circ_diameter)
    M = circles_along_w * circles_along_h

    print(f"Circles that fit along the 140cm width: floor({mat_w} / {circ_diameter}) = {circles_along_w}")
    print(f"Circles that fit along the 110cm height: floor({mat_h} / {circ_diameter}) = {circles_along_h}")
    print(f"Maximum number of circles (M) = {circles_along_w} * {circles_along_h} = {M}")
    print("-" * 40)

    # --- Step 2: Calculate the number of squares (N) in the remaining area ---
    print("Step 2: Calculate the number of squares (N) from the remaining material")
    # The circles are packed into a solid block to maximize leftover rectangular areas.
    packed_circles_w = circles_along_w * circ_diameter
    packed_circles_h = circles_along_h * circ_diameter
    print(f"The {M} circles are packed into a {packed_circles_w}x{packed_circles_h}cm block.")

    # This leaves an L-shaped area, which we divide into two rectangles.
    # Rectangle A: The area along the full width of the material.
    rem_area_A_w = mat_w
    rem_area_A_h = mat_h - packed_circles_h
    sq_in_area_A = math.floor(rem_area_A_w / sq_size) * math.floor(rem_area_A_h / sq_size)

    # Rectangle B: The remaining area next to the packed circles block.
    rem_area_B_w = mat_w - packed_circles_w
    rem_area_B_h = packed_circles_h
    sq_in_area_B = math.floor(rem_area_B_w / sq_size) * math.floor(rem_area_B_h / sq_size)

    N = sq_in_area_A + sq_in_area_B
    print(f"The remaining space is split into two rectangles:")
    print(f" - Area A ({rem_area_A_w}x{rem_area_A_h}cm) fits {sq_in_area_A} squares.")
    print(f" - Area B ({rem_area_B_w}x{rem_area_B_h}cm) fits {sq_in_area_B} squares.")
    print(f"Maximum number of squares (N) = {sq_in_area_A} + {sq_in_area_B} = {N}")
    print("-" * 40)

    # --- Step 3: Calculate the total characters (K) and show the final result ---
    print("Step 3: Calculate the total characters (K) and show the final result")
    k_from_squares = chars_per_sq * N
    k_from_circles = chars_per_circ * M
    K = k_from_squares + k_from_circles

    print("The final calculation for the total number of characters is:")
    print(f"K = (Squares × Chars/Square) + (Circles × Chars/Circle)")
    print(f"K = ({N} × {chars_per_sq}) + ({M} × {chars_per_circ})")
    print(f"K = {k_from_squares} + {k_from_circles}")
    print(f"K = {K}")
    print("-" * 40)

    # Final Answer in N:M:K format
    final_answer = f"{N}:{M}:{K}"
    print(f"The optimal result in the format N:M:K is: {final_answer}")

if __name__ == "__main__":
    solve_engraving_problem()
