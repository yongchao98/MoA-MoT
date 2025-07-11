def solve_poker_strategy():
    """
    Analyzes the poker scenario to determine the optimal strategy for holding AA.
    """
    # 1. Define the game state
    pot = 10.0
    stack = 1000.0
    
    # On the 2s2d2c board, your AA (Aces full of Deuces) beats Villain's KK (Kings full of Deuces).
    # You effectively have the nuts. There are no more cards to change the outcome.

    # 2. Analyze the EV of checking with AA (The Passive Line)
    # If you check with AA, you must do so as part of a wider, unexploitable checking strategy.
    # The simplest GTO approach is to check with your entire range.
    # Villain holds KK and knows your range is {AA, QQ, ..., 33} (1 combo he loses to, 10 he beats).
    # Villain must decide whether to check or bet.
    #   - EV_villain(Check): Wins the pot vs 10 hands, loses vs 1. EV = (10/11) * pot = $9.09
    #   - EV_villain(Bet B): You would only call with AA. Villain wins the pot when you fold (10/11) and loses his bet when you call (1/11).
    #     EV = (10/11)*pot - (1/11)*B = (100 - B) / 11.
    # Comparing these, Villain's EV is always higher if he checks back (for any B>0).
    # Therefore, Villain's best response to you checking your range is to also check.
    # When Villain checks back, you win the pot at showdown.
    ev_check_line = pot

    # 3. Analyze the EV of betting with AA (The Aggressive Line)
    # To get Villain to call a bet, your betting range must contain both value hands (AA) and bluffs.
    # A GTO strategy chooses a bluffing frequency that makes Villain indifferent to calling.
    # Villain's indifference is met when: Bluff_Combos / Value_Combos = Bet / (Pot + Bet)
    # To maximize your profit with AA, you should choose the largest possible bet size, which is all-in.
    bet_size = stack
    
    # Let's confirm you can build a balanced range with this bet size.
    # Required Bluff/Value Ratio = 1000 / (10 + 1000) = 0.99
    # You have 1 value combo (AA) and 10 available bluff combos (QQ-33). You only need ~1 bluff combo
    # to balance your 1 value combo, which is easily done by mixing.
    
    # Since Villain is indifferent, he will mix between calling and folding. We can assume a 50/50 split.
    # EV(Betting AA) = (Prob_Call * Win_if_Called) + (Prob_Fold * Win_if_Folded)
    win_if_called = pot + stack
    win_if_folded = pot
    ev_bet_line = (0.5 * win_if_called) + (0.5 * win_if_folded)
    
    # 4. Compare EVs and Determine the Final Strategy
    # Since EV(Bet) is much greater than EV(Check), you should always bet with AA.
    # The frequency for this action is 100%.
    final_action = "BET"
    # Round sizing and percentage to the nearest even number as per instructions.
    final_sizing = int(round(bet_size / 2.0) * 2.0)
    final_frequency = int(round(100.0 / 2.0) * 2.0)

    print("Step-by-step EV Calculation:")
    print(f"The pot is {pot:.0f} and the effective stack is {stack:.0f}.")
    print("Your AA is the effective nuts against Villain's KK.")
    print(f"EV of a checking strategy with AA: Villain checks back, you win the pot = ${ev_check_line:.0f}")
    print(f"EV of an all-in betting strategy with AA: 0.5 * (Pot + Stack) + 0.5 * Pot")
    print(f"Calculation: 0.5 * ({pot:.0f} + {stack:.0f}) + 0.5 * {pot:.0f} = ${ev_bet_line:.0f}")
    print("\nConclusion: The betting line is vastly more profitable. Therefore, you should always bet all-in with AA.")

    print("\n---")
    print("Final Strategy for AA:")
    print(f"{final_action} {final_sizing} {final_frequency}%")


solve_poker_strategy()