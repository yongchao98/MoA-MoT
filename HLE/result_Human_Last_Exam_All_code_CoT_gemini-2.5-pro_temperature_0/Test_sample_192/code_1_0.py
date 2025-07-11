def solve_poker_problem():
    """
    This function prints the solution to the poker problem based on the GTO analysis.
    The logic is derived from game theory, as explained in the steps above.
    """

    # A) Hero's action with AA
    # Optimal strategy is to bet an amount > $10. We choose the all-in amount.
    # The frequency is 100% as this is a pure strategy.
    action_AA = "BET"
    sizing_AA = 1000
    frequency_AA = 100

    # B) Hero's action with QQ
    # To make the strategy work, we must bet with our entire range, including bluffs.
    action_QQ = "BET"
    sizing_QQ = 1000
    frequency_QQ = 100

    # C) Villain's response
    # Since our bet size ($1000) is greater than $10, and our betting range is 50% value / 50% bluffs,
    # the villain's pot odds are worse than his equity. A perfect player must fold.
    action_villain = "FOLD"
    frequency_villain = 100

    # Print the final answer in the required format.
    # The equation for villain's decision to fold is: P(Hero has bluff) < Bet / (Pot + Bet)
    # 0.5 < 1000 / (10 + 1000)
    # 0.5 < 1000 / 1010
    # 0.5 < 0.99, which is true, so villain folds.
    print(f"A) {action_AA} ${sizing_AA} {frequency_AA}%")
    print(f"B) {action_QQ} ${sizing_QQ} {frequency_QQ}%")
    print(f"C) {action_villain} {frequency_villain}%")

solve_poker_problem()
# <<<A) BET $1000 100% B) BET $1000 100% C) FOLD 100%>>>