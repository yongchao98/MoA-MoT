def solve_rope_puzzle():
    """
    This function solves the magical rope puzzle by demonstrating the
    sequence of cuts that results in a 15cm piece.
    """
    initial_rope = 60
    target_piece = 15
    shrink_factor = 0.75

    # The length of the longer portion after the first cut. This is our answer.
    longer_portion_cut1 = 40

    # Step 1: Perform the first cut.
    shorter_portion_cut1 = initial_rope - longer_portion_cut1
    print(f"Step 1: The initial {initial_rope}cm rope is cut into a longer portion of {longer_portion_cut1}cm and a shorter portion of {shorter_portion_cut1}cm.")
    print(f"{initial_rope} = {longer_portion_cut1} + {shorter_portion_cut1}")
    print("-" * 20)

    # Step 2: The longer portion shrinks.
    shrunken_longer_portion = int(longer_portion_cut1 * shrink_factor)
    print(f"Step 2: The longer portion ({longer_portion_cut1}cm) shrinks by 25%.")
    print(f"{longer_portion_cut1} * {shrink_factor} = {shrunken_longer_portion}")
    print(f"The pieces are now {shorter_portion_cut1}cm and {shrunken_longer_portion}cm long.")
    print("-" * 20)

    # Step 3: Perform the second cut on the 30cm piece.
    piece_for_cut2 = shrunken_longer_portion
    final_piece_1 = piece_for_cut2 / 2
    final_piece_2 = piece_for_cut2 / 2
    print(f"Step 3: The {piece_for_cut2}cm piece is cut into two equal halves.")
    print(f"Since there is no 'longer portion', no shrinkage occurs.")
    print(f"{piece_for_cut2} / 2 = {int(final_piece_1)}")
    print("-" * 20)

    # Conclusion
    print(f"Result: A piece of {int(final_piece_1)}cm is obtained.")
    print(f"\nThe required length of the longer portion after the first cut is {longer_portion_cut1}cm.")

solve_rope_puzzle()