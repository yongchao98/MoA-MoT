def solve_poker_gto_question():
    """
    Analyzes the reasons for betting in poker under the assumption of
    Game-Theoretically Optimal (GTO) play.
    """
    reasons = {
        3: "Denying equity to drawing hands.",
        4: "Gaining information about your opponent's hand.",
        5: "Avoiding revealing your own hand in a showdown."
    }

    disappearing_reasons = []
    explanation = []

    # Analysis of reason 3
    # Denying equity is a core part of GTO. A balanced betting range forces draws to fold or pay incorrectly.
    # This reason does NOT disappear.

    # Analysis of reason 4
    # Against a GTO player, their entire strategy is known. A bet doesn't gain new, exploitable
    # information because their response (call/fold/raise) is part of a known, balanced probability.
    disappearing_reasons.append(4)
    explanation.append("Reason 4 (Gaining Information) disappears because a GTO opponent's strategy is already known. Their reaction to a bet does not provide new exploitable information, but is a pre-defined probabilistic response.")

    # Analysis of reason 5
    # A GTO player already knows your complete strategy. Showing one hand at showdown
    # doesn't reveal anything they can exploit, as a GTO strategy is unexploitable by definition.
    disappearing_reasons.append(5)
    explanation.append("Reason 5 (Avoiding Revealing Your Hand) disappears because your GTO strategy is also known and is unexploitable. There is no strategic information to hide from an opponent who cannot exploit it.")

    print("In a game where both players use a game-theoretically optimal strategy, the following subtle reasons for betting disappear:")
    print("-" * 80)
    for text in explanation:
        print(text)
    print("-" * 80)
    
    # Format the final numeric answer as a string
    final_answer_numbers = " and ".join(map(str, sorted(disappearing_reasons)))
    print(f"Therefore, reasons {final_answer_numbers} disappear.")

solve_poker_gto_question()