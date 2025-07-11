def solve_pll_stickers():
    """
    This function explains and calculates the minimum number of non-top-facing stickers
    required to identify a PLL case on a 3x3 Rubik's Cube.
    """

    # There are 21 distinct PLL cases to differentiate.
    num_pll_cases = 21

    # The analysis focuses on the "worst-case scenario" - distinguishing the most
    # similar PLL cases, which are typically found within the four G-perms.
    
    # Let's consider distinguishing two similar G-perms (e.g., Ga and Gd).
    # Both can be oriented to show "headlights" on one face (e.g., the right face).
    # "Headlights" are when the side stickers of the two corners on a face match.
    # To confirm you have headlights, you must observe these two stickers.
    stickers_for_headlights = 2
    
    # After observing the headlights, the two G-perms can still appear identical.
    # The feature that distinguishes them is the location of a 1x2 "block"
    # (an adjacent edge and corner whose side stickers match).
    # To find this block, you need to look at a different face. For example, to check
    # for a block on the front-left, you need to see both the corner and edge stickers there.
    stickers_for_block = 2
    
    # Therefore, to distinguish these worst-case scenarios, one must see the
    # stickers for the headlights AND the stickers for the block.
    min_stickers_required = stickers_for_headlights + stickers_for_block
    
    print("To solve this puzzle, we analyze the worst-case scenario for PLL recognition.")
    print(f"There are a total of {num_pll_cases} PLL cases to identify.")
    print("The hardest cases to distinguish are often G-permutations.")
    print("\nTo differentiate between two very similar G-perms:")
    print(f"1. You first identify a pattern like 'headlights', which requires observing at least {stickers_for_headlights} stickers.")
    print(f"2. To distinguish between them, you must then locate a 'block' on another face, which requires observing another {stickers_for_block} stickers.")
    
    # Final equation and result
    print(f"\nThe minimum number of stickers is the sum of these observations:")
    print(f"{stickers_for_headlights} (for headlights) + {stickers_for_block} (for the block) = {min_stickers_required}")
    print(f"\nThus, you must see a minimum of {min_stickers_required} non-top-facing stickers to fully identify any PLL case.")

solve_pll_stickers()