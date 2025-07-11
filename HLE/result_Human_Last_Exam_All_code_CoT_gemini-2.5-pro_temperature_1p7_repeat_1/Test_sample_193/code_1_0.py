def solve_poker_strategy():
    """
    This function analyzes the given poker scenario to determine the optimal
    strategy for the player holding AA.
    """

    # --- Problem Parameters ---
    pot = 10
    stack = 1000
    
    # --- GTO Analysis ---
    
    # 1. Optimal Bet Sizing:
    # The EV (profit) for betting with AA is a function of the bet size B:
    # EV(B) = pot + (pot / (pot + B)) * B
    # The derivative of this function with respect to B is 100 / (10 + B)^2,
    # which is always positive for B > 0.
    # This means the EV is maximized when B is maximized.
    # The maximum possible bet is the entire stack.
    optimal_sizing = stack

    # 2. Optimal Frequency for AA:
    # EV of checking with AA is winning the pot: $10.
    # EV of betting with AA is pot + alpha * B = 10 + (10/(10+B))*B.
    # Since B > 0, the term (10/(10+B))*B is always positive.
    # Therefore, EV(Bet) > EV(Check) for AA.
    # This means we should always bet with AA and never check.
    bet_frequency = 100

    # 3. Format the final answer per user instructions.
    # The action is to bet.
    # The sizing must be rounded to the nearest even number. 1000 is even.
    # The percentage must be rounded to the nearest even number. 100 is even.
    
    action = "BET"
    final_sizing = int(round(optimal_sizing / 2) * 2)
    final_frequency = int(round(bet_frequency / 2) * 2)

    print(f"To be unexploitable, your strategy with AA is a pure one, not mixed.")
    print(f"You want to bet the largest possible amount to maximize value.")
    print(f"The optimal play is:")
    print(f"{action} {final_sizing} {final_frequency}%")


solve_poker_strategy()