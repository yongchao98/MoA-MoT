def solve_pll_sticker_count():
    """
    Calculates the minimum number of non-top-facing stickers needed
    to identify any PLL case on a 3x3 Rubik's Cube.
    """
    # The number of stickers to identify a 1x2 solved "block"
    # (one corner side-sticker and one edge side-sticker).
    stickers_for_block = 2

    # The number of stickers to identify "headlights" on an adjacent face
    # (the two side-stickers of the corner pieces).
    stickers_for_headlights = 2

    # The number of stickers for the "tie-breaker" piece.
    # This is the edge piece located between the headlights. Its state
    # distinguishes the T-perm from the four G-perms.
    stickers_for_tiebreaker_edge = 1

    # The total is the sum of stickers needed to see the ambiguous pattern
    # plus the one sticker needed to resolve the ambiguity.
    total_stickers_needed = stickers_for_block + stickers_for_headlights + stickers_for_tiebreaker_edge

    print("On a 3x3 Rubik's Cube, identifying the PLL case requires looking at the side stickers of the top layer.")
    print("The minimum number of stickers one must see is determined by the most difficult cases to distinguish.")
    print("\nThe hardest cases to tell apart are the T-perm and the four G-perms.")
    print("To guarantee a correct identification, a solver must observe the following:")
    print(f"\n1. The two stickers forming a 1x2 'block'.")
    print(f"2. The two stickers forming 'headlights' on an adjacent face.")
    print(f"3. The single sticker on the edge piece between the headlights, which acts as the tie-breaker.")
    
    print("\nTherefore, the minimum number of stickers required is the sum of stickers for these features:")
    print(f"{stickers_for_block} + {stickers_for_headlights} + {stickers_for_tiebreaker_edge} = {total_stickers_needed}")

solve_pll_sticker_count()