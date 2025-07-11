def solve_packing_problem():
    """
    Finds the optimal number of squares (N) and circles (M) to maximize
    the total number of engraved characters (K) from a single sheet.
    """
    
    # Define problem constants
    sheet_dims = [(140, 110), (110, 140)] # Two possible orientations
    square_dim = 10
    circle_bbox = 40
    chars_per_square = 4
    chars_per_circle = 999
    
    max_k = -1
    best_config = {}

    # Iterate through both sheet orientations
    for W, H in sheet_dims:
        
        # Iterate through possible numbers of circle columns
        for m_w in range(W // circle_bbox + 1):
            
            # Iterate through possible numbers of circle rows
            for m_h in range(H // circle_bbox + 1):
                
                # Number of circles for this layout
                M = m_w * m_h
                
                # Dimensions of the rectangular block reserved for circles
                circles_block_w = m_w * circle_bbox
                circles_block_h = m_h * circle_bbox
                
                # The remaining L-shaped area can be perfectly tiled by squares because
                # all dimensions are multiples of the square's dimension.
                # Decompose the L-shape into two rectangles:
                
                # 1. Rectangle to the right of the circle block
                rem1_w = W - circles_block_w
                rem1_h = H # Use full height H, as the area below this is handled by rem2
                
                # 2. Rectangle below the circle block (and not overlapping with rem1)
                rem2_w = circles_block_w
                rem2_h = H - circles_block_h

                # Calculate number of squares in each remaining rectangle
                n1 = (rem1_w // square_dim) * (rem1_h // square_dim)
                n2 = (rem2_w // square_dim) * (rem2_h // square_dim)
                
                N = n1 + n2

                # Calculate total characters
                K = chars_per_square * N + chars_per_circle * M
                
                if K > max_k:
                    max_k = K
                    best_config = {'N': N, 'M': M, 'K': K}
                    
    # Final results
    N = best_config['N']
    M = best_config['M']
    K = best_config['K']
    
    print(f"The optimal combination is to produce {N} squares and {M} circles.")
    print("\nEach square has 4 characters, and each circle has 999 characters.")
    print("The final calculation is:")
    print(f"K = ({N} * {chars_per_square}) + ({M} * {chars_per_circle})")
    print(f"K = {N * chars_per_square} + {M * chars_per_circle}")
    print(f"K = {K}")

    print(f"\nThe final answer in the format N:M:K is: {N}:{M}:{K}")

solve_packing_problem()
<<<58:6:6226>>>