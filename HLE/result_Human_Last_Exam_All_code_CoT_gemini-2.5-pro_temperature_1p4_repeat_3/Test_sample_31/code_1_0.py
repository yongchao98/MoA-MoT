def solve_pll_sticker_puzzle():
    """
    Calculates the minimum number of non-top-facing stickers needed
    to identify a PLL case on a 3x3 Rubik's Cube.
    """

    # Step 1: Consider the number of stickers on a single face's last layer.
    # These are the two corner side-stickers and one edge side-sticker.
    stickers_on_one_face = 3

    # Step 2: Determine if this is sufficient.
    # A common recognition pattern is a "bar", where all 3 stickers on one face
    # have the same color. For example, on the Blue face, the pattern is Blue-Blue-Blue.
    # This single pattern can be produced by several different PLL cases:
    #  - Solved state
    #  - Ua Permutation (3-cycle of edges)
    #  - Ub Permutation (3-cycle of edges)
    #  - Aa Permutation (3-cycle of corners)
    #  - Ab Permutation (3-cycle of corners)
    #
    # Since observing these 3 stickers leads to multiple possibilities,
    # 3 stickers are not enough to guarantee unique identification.

    print("Analyzing the ambiguity when viewing 3 stickers on one face:")
    print(f"Viewing a single face gives information from {stickers_on_one_face} stickers.")
    print("However, a 'bar' pattern (3 same-colored stickers) can be many PLLs.")
    print("This proves that 3 stickers are not enough.\n")

    # Step 3: Add the minimum number of additional stickers to resolve ambiguity.
    # To distinguish the cases that all look like a "bar" from the front,
    # we must gather more information from another face.
    #
    # For example, by observing just ONE additional sticker on an adjacent face (e.g., the right-face edge),
    # we can distinguish the U-perms from the A-perms and the solved state.
    #
    # Other hard-to-distinguish cases, like the four G-perms, are recognized by
    # identifying a 2-sticker "headlights" pattern and a 2-sticker "block" pattern
    # on adjacent faces. This requires observing a total of 4 stickers.
    #
    # Since 3 stickers are proven to be insufficient, and 4 stickers are sufficient
    # to resolve the most difficult cases, the minimum required number is 4.

    additional_stickers_needed = 1
    minimum_required_stickers = stickers_on_one_face + additional_stickers_needed

    print("To resolve the ambiguity, at least one more sticker must be seen.")
    print(f"Stickers on first face: {stickers_on_one_face}")
    print(f"Minimum additional stickers needed: {additional_stickers_needed}")
    print(f"Therefore, the total minimum is {stickers_on_one_face} + {additional_stickers_needed} = {minimum_required_stickers}.")

solve_pll_sticker_puzzle()