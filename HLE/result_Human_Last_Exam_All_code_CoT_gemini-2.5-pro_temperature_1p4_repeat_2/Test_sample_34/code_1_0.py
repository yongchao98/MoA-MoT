def solve_poker_question():
    """
    This function determines the worst suited jack to open from the button
    in a 100BB rake-free cash game and prints the result.

    The ranking of suited jacks (JXs) from best to worst is based on the kicker card.
    Higher kickers are better due to less domination and better straight potential.
    - Best: AJs, KJs, QJs, JTs (strong broadway cards and connectors)
    - Middle: J9s, J8s, J7s (suited connectors and one-gappers)
    - Weakest: J6s, J5s, J4s, J3s, J2s (low kickers, poor connectivity)

    In a standard GTO button opening range, all suited jacks are typically included
    as a raise. The one with the lowest value, and therefore the "worst" one you should
    still open, is the one with the lowest kicker.
    """
    worst_suited_jack = "J2s"
    print(worst_suited_jack)

solve_poker_question()