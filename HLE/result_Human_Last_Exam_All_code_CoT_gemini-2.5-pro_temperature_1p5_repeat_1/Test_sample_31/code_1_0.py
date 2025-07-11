def solve_pll_sticker_puzzle():
    """
    Calculates and explains the minimum number of non-top-facing stickers
    that must be seen to fully identify the PLL case.
    """

    # Step 1: Determine the Edge Permutation (Pe)
    # The last layer has 4 edge pieces. A permutation of 4 items is fully
    # determined if we know the destination of 3 of them.
    # The destination of an edge piece can be determined by observing the color
    # of its single non-top-facing sticker.
    # For example, by looking at the sticker on the front face of the front-edge position,
    # we know which of the 4 edge pieces has moved there.
    stickers_for_edges = 3
    print(f"Stickers to determine the full edge permutation: {stickers_for_edges}")

    # Step 2: Determine the Corner Permutation (Pc)
    # The last layer has 4 corner pieces. Each corner has 2 non-top-facing stickers.
    # To uniquely identify a corner piece, we must see both of its side stickers.

    # We can determine the location of the first two corner pieces.
    stickers_for_corner_1 = 2
    print(f"Stickers to identify the piece in the first corner position: {stickers_for_corner_1}")
    stickers_for_corner_2 = 2
    print(f"Stickers to identify the piece in the second corner position: {stickers_for_corner_2}")

    # Step 3: Use the Parity rule to find the rest.
    # Once we know where 2 of the 4 corners are, the remaining 2 corners must be
    # in the remaining 2 positions. There are only two possibilities for this:
    #   1. Corner 3 is in position 3, and corner 4 is in position 4.
    #   2. Corner 3 is in position 4, and corner 4 is in position 3 (a swap).
    #
    # These two possibilities for the full corner permutation have opposite parity.
    # Crucially, the parity of the corner permutation MUST match the parity of the
    # edge permutation, which we already determined in Step 1.
    # Therefore, we can deduce which of the two possibilities is the correct one
    # without looking at any more stickers.
    stickers_for_remaining_corners = 0
    print(f"Stickers for the remaining two corners (using parity rule): {stickers_for_remaining_corners}")
    
    # Step 4: Calculate the total minimum number of stickers.
    total_stickers = stickers_for_edges + stickers_for_corner_1 + stickers_for_corner_2 + stickers_for_remaining_corners
    
    print("\n---")
    print("The final calculation is:")
    print(f"Total Stickers = {stickers_for_edges} (edges) + {stickers_for_corner_1} (corner 1) + {stickers_for_corner_2} (corner 2) + {stickers_for_remaining_corners} (parity) = {total_stickers}")

solve_pll_sticker_puzzle()