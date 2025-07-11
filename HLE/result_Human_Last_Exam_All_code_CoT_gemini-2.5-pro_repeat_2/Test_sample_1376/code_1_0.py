import math

def solve_tomb_engraving():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # Define constants based on the problem description
    MATERIAL_W = 140
    MATERIAL_H = 110
    SQUARE_SIDE = 10
    CIRCLE_DIAMETER = 40
    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    best_config = {'N': 0, 'M': 0, 'K': 0, 'description': 'Initial'}

    # --- Strategy 1: Only Squares ---
    n1 = (MATERIAL_W // SQUARE_SIDE) * (MATERIAL_H // SQUARE_SIDE)
    m1 = 0
    k1 = n1 * CHARS_PER_SQUARE + m1 * CHARS_PER_CIRCLE
    if k1 > best_config['K']:
        best_config = {'N': n1, 'M': m1, 'K': k1, 'description': 'Only Squares'}

    # --- Strategy 2: Grid-Packed Circles + Squares in Remainder ---
    # Circles are arranged in a 3x2 grid
    m2 = (MATERIAL_W // CIRCLE_DIAMETER) * (MATERIAL_H // CIRCLE_DIAMETER) # 3 * 2 = 6 circles
    circles_w_used = (MATERIAL_W // CIRCLE_DIAMETER) * CIRCLE_DIAMETER # 3 * 40 = 120
    circles_h_used = (MATERIAL_H // CIRCLE_DIAMETER) * CIRCLE_DIAMETER # 2 * 40 = 80
    
    # The remaining area is an L-shape, which we can split into two rectangles
    rem_rect1_w = MATERIAL_W - circles_w_used # 140 - 120 = 20
    rem_rect1_h = MATERIAL_H # 110
    rem_rect2_w = circles_w_used # 120
    rem_rect2_h = MATERIAL_H - circles_h_used # 110 - 80 = 30
    
    squares_in_rect1 = (rem_rect1_w // SQUARE_SIDE) * (rem_rect1_h // SQUARE_SIDE) # (20/10)*(110/10) = 2*11=22
    squares_in_rect2 = (rem_rect2_w // SQUARE_SIDE) * (rem_rect2_h // SQUARE_SIDE) # (120/10)*(30/10) = 12*3=36
    n2 = squares_in_rect1 + squares_in_rect2 # 22 + 36 = 58
    k2 = n2 * CHARS_PER_SQUARE + m2 * CHARS_PER_CIRCLE
    if k2 > best_config['K']:
        best_config = {'N': n2, 'M': m2, 'K': k2, 'description': 'Grid-Packed Circles'}

    # --- Strategy 3: Densely Packed Circles (Max M) ---
    # A 3+3+3 hexagonal packing of 9 circles fits in a 140 x 109.28cm box.
    # This leaves no space for any 10x10cm squares.
    m3 = 9
    n3 = 0
    k3 = n3 * CHARS_PER_SQUARE + m3 * CHARS_PER_CIRCLE
    if k3 > best_config['K']:
        best_config = {'N': n3, 'M': m3, 'K': k3, 'description': 'Max Dense-Packed Circles'}

    # --- Strategy 4: Densely Packed Circles (Sub-Max M) ---
    # A 3+2+3 hexagonal packing of 8 circles fits in a 120 x 109.28cm box.
    # This leaves a 20 x 109.28cm rectangular area for squares.
    m4 = 8
    rem_w = MATERIAL_W - 120 # 140 - 120 = 20
    rem_h = 109.28 # Height of the packed circles area
    n4 = (rem_w // SQUARE_SIDE) * (math.floor(rem_h) // SQUARE_SIDE) # (20/10)*(109/10) = 2*10=20
    k4 = n4 * CHARS_PER_SQUARE + m4 * CHARS_PER_CIRCLE
    if k4 > best_config['K']:
        best_config = {'N': n4, 'M': m4, 'K': k4, 'description': 'Sub-Max Dense-Packed Circles'}

    # --- Final Result ---
    final_N = best_config['N']
    final_M = best_config['M']
    final_K = best_config['K']

    print(f"The optimal configuration is to produce {final_N} squares and {final_M} circles.")
    print("This maximizes the total number of engraved characters.")
    print("\nFinal Calculation:")
    print(f"{final_N} squares * {CHARS_PER_SQUARE} chars/square + {final_M} circles * {CHARS_PER_CIRCLE} chars/circle = {final_K} total characters")
    
    print("\nAnswer in N:M:K format:")
    print(f"{final_N}:{final_M}:{final_K}")


solve_tomb_engraving()