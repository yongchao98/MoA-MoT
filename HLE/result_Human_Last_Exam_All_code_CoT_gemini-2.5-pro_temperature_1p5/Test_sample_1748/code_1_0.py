def solve_rope_puzzle():
    """
    This function solves the magical rope puzzle by working backwards
    and printing the steps of the solution.
    """

    # Start with the final desired piece and work backwards.
    final_piece = 15  # cm

    # Step 1: Analyze the second cut.
    # To get a 15cm piece, we must cut a rope of length L2.
    # To prevent the 15cm piece from shrinking, it must be the shorter or equal part.
    # So, 15 <= L2 - 15, which means L2 >= 30.
    # The simplest case is cutting a 30cm rope in half.
    rope_for_cut2 = 30  # cm

    # Step 2: Analyze the first cut.
    # The 30cm rope for the second cut must come from the first cut on the 60cm rope.
    # Let the longer portion of the first cut be x. The shorter portion is 60-x.
    # The longer portion shrinks by 25% (multiplied by 0.75).
    # We need the shrunken longer portion to be 30cm.
    # So, x * 0.75 = 30
    shrink_factor = 0.75
    longer_portion_cut1 = rope_for_cut2 / shrink_factor

    # Step 3: Present the solution step-by-step.
    initial_rope = 60
    shorter_portion_cut1 = initial_rope - longer_portion_cut1

    print("Step-by-step solution:")
    print(f"1. The first cut divides the {initial_rope}cm rope into a longer portion and a shorter portion.")
    print(f"   The longer portion is {int(longer_portion_cut1)}cm and the shorter portion is {int(shorter_portion_cut1)}cm.")
    print(f"   Equation: {initial_rope}cm -> {int(longer_portion_cut1)}cm + {int(shorter_portion_cut1)}cm")
    print("")
    print("2. The longer portion shrinks by 25%.")
    print(f"   Its new length is {int(longer_portion_cut1)} * {shrink_factor} = {int(rope_for_cut2)}cm.")
    print(f"   Equation: {int(longer_portion_cut1)} * (1 - 0.25) = {int(rope_for_cut2)}")
    print("")
    print("3. This 30cm piece is then cut in half for the second cut.")
    print(f"   This produces two {final_piece}cm pieces. Since they are equal length, no shrinking occurs.")
    print(f"   Equation: {int(rope_for_cut2)}cm -> {final_piece}cm + {final_piece}cm")
    print("")
    print("The question asks for the length of the longer portion after the first cut.")
    print(f"This is the length before shrinking.")
    print(f"The final answer is: {int(longer_portion_cut1)}")


solve_rope_puzzle()
<<<40>>>