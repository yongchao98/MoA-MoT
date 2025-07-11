def solve_poker_problem():
    """
    This function solves the poker GTO problem based on the provided scenario.

    The logic is as follows:
    1.  Analyze Hand Strength: The board is 2s2d2c3h3c (Threes full of Twos).
        Neither Hero's range (AA, QQ) nor Villain's known hand (KK) can improve
        on the board. Therefore, any showdown is a guaranteed chop. All parties know this.

    2.  Analyze Villain's Response to a Bet:
        - Pot = $10.
        - If Hero bets B and Villain folds, Villain's EV = $0.
        - If Hero bets B and Villain calls, the pot is $10 + 2B. They chop it.
          Villain receives ($10 + 2B) / 2 = $5 + B.
          Villain's profit = ($5 + B) - B = $5.
        - Since $5 > $0, Villain will always call any bet.

    3.  Analyze Hero's Action:
        - If Hero checks, they chop the $10 pot. Hero's profit = $0.
        - If Hero bets B, Villain calls. They chop the pot.
          Hero's profit = ($5 + B) - B = $5.
        - Since $5 > $0, Hero should always bet.

    4.  Determine Bet Sizing:
        - The profit is $5 for any bet size B > 0.
        - A logical strategy minimizes risk for the same reward.
        - We choose a small, non-zero bet size. The problem asks to round to an
          even number, so $2 is a perfect choice.

    5.  Construct the final answer string.
    """

    # Case A: Hero holds AA
    # Strategy is to bet 100% of the time.
    bet_sizing_a = 2
    action_a = f"BET ${bet_sizing_a} 100%"

    # Case B: Hero holds QQ
    # The situation is identical to holding AA.
    bet_sizing_b = 2
    action_b = f"BET ${bet_sizing_b} 100%"

    # Case C: Villain's response
    # As determined above, Villain's only profitable move is to call.
    action_c = "CALL 100%"

    # Format the final output string
    final_answer = f"A) {action_a} B) {action_b} C) {action_c}"
    print(final_answer)


solve_poker_problem()