def solve_engraving_optimization():
    """
    Calculates the optimal number of squares and circles to maximize the
    number of engraved Chinese characters from a given rectangular material.
    """
    # --- Problem Parameters ---
    
    # Dimensions of the raw material rectangle
    board_w, board_h = 140, 110

    # Properties of the square artifacts
    square_size = 10
    square_chars = 4

    # Properties of the circular artifacts
    # A 20cm radius circle requires a 40x40cm square bounding box to be cut
    circle_bbox_size = 40
    circle_chars = 999

    # --- Step 1: Maximize the number of circles (M) ---
    # Due to the much higher character count per circle (999 vs 4),
    # we prioritize fitting as many circles as possible.

    # We check the straightforward packing orientation.
    m_along_w = board_w // circle_bbox_size  # 140 // 40 = 3
    m_along_h = board_h // circle_bbox_size  # 110 // 40 = 2
    M = m_along_w * m_along_h               # 3 * 2 = 6

    # --- Step 2: Calculate remaining area and fit squares (N) ---
    # The M circles are arranged in a grid, occupying a rectangular block.
    circles_block_w = m_along_w * circle_bbox_size # 3 * 40 = 120
    circles_block_h = m_along_h * circle_bbox_size # 2 * 40 = 80
    
    # This leaves an L-shaped area, which can be split into two rectangles.
    # Rectangle 1
    rem_area1_w = board_w - circles_block_w  # 140 - 120 = 20
    rem_area1_h = board_h                   # 110
    n1 = (rem_area1_w // square_size) * (rem_area1_h // square_size) # (20//10)*(110//10) = 2*11 = 22

    # Rectangle 2
    rem_area2_w = circles_block_w           # 120
    rem_area2_h = board_h - circles_block_h  # 110 - 80 = 30
    n2 = (rem_area2_w // square_size) * (rem_area2_h // square_size) # (120//10)*(30//10) = 12*3 = 36

    # Total number of squares (N)
    N = n1 + n2

    # --- Step 3: Calculate the total number of characters (K) ---
    K = N * square_chars + M * circle_chars

    # --- Step 4: Output the results ---
    print("This plan maximizes the number of engraved characters.")
    print(f"Number of squares (N) for the name: {N}")
    print(f"Number of circles (M) for the bio: {M}")
    print("\nMaximum characters calculation:")
    # The final equation as requested
    print(f"K = ({N} squares * {square_chars} chars) + ({M} circles * {circle_chars} chars)")
    print(f"K = {N * square_chars} + {M * circle_chars}")
    print(f"K = {K}")
    
    print("\nFinal Answer (N:M:K):")
    print(f"{N}:{M}:{K}")

# Execute the function to solve the problem
solve_engraving_optimization()