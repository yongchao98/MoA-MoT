def solve_poker_gto_question():
    """
    Analyzes poker betting reasons under the assumption of Game-Theoretically Optimal (GTO) play.
    """

    reasons = {
        3: "Denying equity to drawing hands.",
        4: "Gaining information about your opponent's hand.",
        5: "Avoiding revealing your own hand in a showdown."
    }

    # In GTO, the opponent's strategy is a known, perfectly balanced model.
    # Therefore, betting to "gain information" is a moot point. The information is already priced in.
    # A GTO strategy is unexploitable by definition, even if fully known.
    # Therefore, hiding your hand to prevent being figured out is not a concern.
    # Denying equity, however, is a core mathematical component of constructing GTO betting ranges.

    disappearing_reasons_numbers = [4, 5]

    print("Analyzing the subtler reasons for betting when both players use a Game-Theoretically Optimal (GTO) strategy:")
    print("-" * 80)
    print(f"Reason 3: '{reasons[3]}' -> This reason REMAINS. Denying equity is a fundamental mathematical component of GTO strategy.")
    print(f"Reason {disappearing_reasons_numbers[0]}: '{reasons[4]}' -> This reason DISAPPEARS. In a GTO vs. GTO context, your opponent's strategy is a known quantity (the GTO model itself). You don't need to bet to 'gain' information; their responses are already accounted for in the equilibrium.")
    print(f"Reason {disappearing_reasons_numbers[1]}: '{reasons[5]}' -> This reason DISAPPEARS. A GTO strategy is unexploitable by definition. Revealing a single hand at showdown does not make your overall, perfectly balanced strategy vulnerable.")
    print("-" * 80)
    
    # Final conclusion as per the prompt's instructions
    print(f"The logic shows that subtler reasons {disappearing_reasons_numbers[0]} and {disappearing_reasons_numbers[1]} disappear under GTO assumptions.")


solve_poker_gto_question()