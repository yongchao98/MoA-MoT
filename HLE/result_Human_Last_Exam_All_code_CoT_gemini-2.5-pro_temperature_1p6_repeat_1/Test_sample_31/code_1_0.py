def solve_pll_sticker_puzzle():
    """
    Calculates the minimum number of non-top-facing stickers that must be
    seen to guarantee identification of any PLL case.

    The reasoning is based on a worst-case analysis to find a number N
    such that viewing any N stickers is sufficient to identify the case.
    """

    print("### Step 1: Identifying the Worst-Case for Recognition ###")
    print("To guarantee identification, we must consider the hardest-to-distinguish cases.")
    print("This occurs when we try to differentiate the 'solved' state from a PLL case that moves the fewest pieces, as they will share the most identical sticker positions.")
    print("The U-Permutation moves the fewest pieces: it cycles 3 edges, leaving all 4 corners and the remaining 1 edge stationary.\n")

    # Cube and piece parameters
    total_corners = 4
    total_edges = 4
    side_stickers_per_corner = 2
    side_stickers_per_edge = 1

    # U-Permutation piece movement
    moved_edges = 3

    print("### Step 2: Counting Stickers on Stationary Pieces for a U-Perm ###")
    # Calculate how many pieces are stationary
    stationary_corners = total_corners
    stationary_edges = total_edges - moved_edges
    print(f"A U-Perm leaves {stationary_corners} corners and {stationary_edges} edge stationary.")

    # Calculate the number of stickers on those stationary pieces
    stickers_on_stationary_corners = stationary_corners * side_stickers_per_corner
    stickers_on_stationary_edges = stationary_edges * side_stickers_per_edge
    total_stationary_stickers = stickers_on_stationary_corners + stickers_on_stationary_edges
    print(f"The number of non-top stickers on the {stationary_corners} stationary corners is: {stationary_corners} * {side_stickers_per_corner} = {stickers_on_stationary_corners}.")
    print(f"The number of non-top stickers on the {stationary_edges} stationary edge is: {stationary_edges} * {side_stickers_per_edge} = {stickers_on_stationary_edges}.")
    print(f"Total number of stickers that remain unchanged in a U-Perm case is: {stickers_on_stationary_corners} + {stickers_on_stationary_edges} = {total_stationary_stickers}.\n")

    print("### Step 3: Applying the Pigeonhole Principle ###")
    print(f"If a solver looks at only {total_stationary_stickers} stickers, it is possible they chose only stickers that are unaffected by the U-Perm.")
    print("In that scenario, the cube would appear identical to the solved state, and the PLL case could not be identified.")
    print("To guarantee seeing at least one sticker on a moved piece, a solver must look at one more sticker than the total number of stationary stickers.\n")

    # The final calculation
    min_stickers_for_guaranteed_id = total_stationary_stickers + 1
    print("### Conclusion ###")
    print(f"The minimum number of stickers that must be seen is {total_stationary_stickers} + 1 = {min_stickers_for_guaranteed_id}.")


solve_pll_sticker_puzzle()
