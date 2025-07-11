import math

def solve_engraving_maximization():
    """
    Solves the problem of maximizing the number of engraved characters
    by optimally cutting squares and circles from a rectangle.
    """
    # Define the dimensions and character counts
    RECT_W = 140
    RECT_H = 110
    SQUARE_SIDE = 10
    CIRCLE_DIAMETER = 40  # From 20cm radius
    NAME_CHARS_PER_SQUARE = 4
    BIO_CHARS_PER_CIRCLE = 999

    print("Step 1: Determine the strategy.")
    print(f"To maximize K = {NAME_CHARS_PER_SQUARE}*N + {BIO_CHARS_PER_CIRCLE}*M, we prioritize M due to its high coefficient.\n")

    # --- Calculation for one orientation (140x110) ---
    # We will verify with the other orientation (110x140) if necessary,
    # but the result is the same.
    
    print("Step 2: Calculate the maximum number of circles (M).")
    # Pack as many 40x40cm bounding boxes for circles as possible.
    cols_M = RECT_W // CIRCLE_DIAMETER
    rows_M = RECT_H // CIRCLE_DIAMETER
    M = cols_M * rows_M
    print(f"The {RECT_W}x{RECT_H}cm rectangle can fit {cols_M} circles along its width and {rows_M} along its height.")
    print(f"Max number of circles (M) = {cols_M} * {rows_M} = {M}\n")

    # Area occupied by the circles' bounding boxes
    used_W = cols_M * CIRCLE_DIAMETER
    used_H = rows_M * CIRCLE_DIAMETER
    
    print("Step 3: Calculate the number of squares (N) in the remaining area.")
    print(f"The {M} circles occupy a {used_W}x{used_H}cm area, leaving an L-shaped space.")

    # Decompose the L-shaped remaining area into two rectangles.
    # Method 1 of decomposition:
    rem_area1_W = RECT_W - used_W
    rem_area1_H = RECT_H
    rem_area2_W = used_W
    rem_area2_H = RECT_H - used_H
    
    # Calculate how many squares fit in each rectangular part
    squares_in_area1 = (rem_area1_W // SQUARE_SIDE) * (rem_area1_H // SQUARE_SIDE)
    squares_in_area2 = (rem_area2_W // SQUARE_SIDE) * (rem_area2_H // SQUARE_SIDE)
    
    N = squares_in_area1 + squares_in_area2
    
    print(f"The remaining space is split into a {rem_area1_W}x{rem_area1_H}cm area and a {rem_area2_W}x{rem_area2_H}cm area.")
    print(f"Number of squares in first area: ({rem_area1_W}//{SQUARE_SIDE}) * ({rem_area1_H}//{SQUARE_SIDE}) = {squares_in_area1}")
    print(f"Number of squares in second area: ({rem_area2_W}//{SQUARE_SIDE}) * ({rem_area2_H}//{SQUARE_SIDE}) = {squares_in_area2}")
    print(f"Total number of squares (N) = {squares_in_area1} + {squares_in_area2} = {N}\n")
    
    print("Step 4: Calculate the total number of engraved characters (K).")
    K = NAME_CHARS_PER_SQUARE * N + BIO_CHARS_PER_CIRCLE * M
    
    # To satisfy the prompt's request: "output each number in the final equation"
    k_from_n = NAME_CHARS_PER_SQUARE * N
    k_from_m = BIO_CHARS_PER_CIRCLE * M
    
    print(f"The final equation is: K = ({NAME_CHARS_PER_SQUARE} * N) + ({BIO_CHARS_PER_CIRCLE} * M)")
    print(f"K = ({NAME_CHARS_PER_SQUARE} * {N}) + ({BIO_CHARS_PER_CIRCLE} * {M})")
    print(f"K = {k_from_n} + {k_from_m}")
    print(f"K = {K}\n")
    
    print("Final Answer:")
    # Output the result in the specified format N:M:K
    final_answer = f"{N}:{M}:{K}"
    print(final_answer)
    
solve_engraving_maximization()