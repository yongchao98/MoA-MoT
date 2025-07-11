import math

def solve_poker_strategy():
    """
    Calculates the optimal poker strategy for the given scenario.
    """
    # Game state variables
    pot = 10
    stack = 1000
    
    # Hero's range analysis
    # Hero has AA (value) and 10 pairs from QQ-33 (bluffs/losers).
    # Each has 6 combos, but the ratio is what matters.
    # Value hands (beat KK): AA
    # Losing hands (lose to KK): QQ, JJ, TT, 99, 88, 77, 66, 55, 44, 33
    
    # --- Step 1: Calculate EV of checking with AA ---
    # If Hero checks with a mix of AA and losing hands, Villain (with KK)
    # sees a range that he mostly beats. He should bet for value.
    # To maximize value, Villain should bet an amount that Hero's losing hands
    # will fold to, but that Hero's AA will call/raise.
    # Let's assume Villain makes a bet. If Hero has AA, he can trap and win Villain's bet.
    # Hero's AA is the nuts. If Villain bets, Hero can raise all-in. Villain, knowing
    # a raise only comes from AA, will fold to the raise. Hero wins the pot + Villain's initial bet.
    # Villain's optimal bet size when facing a range of nuts and bluff-catchers he beats
    # is to bet the minimum possible amount, to lose the minimum to the nuts while
    # getting folds from the bluff-catchers.
    # So, Villain min-bets, Hero (AA) raises, Villain folds. Hero wins the pot.
    ev_check_aa = pot
    print("Poker Strategy Calculation:")
    print("---------------------------------")
    print("1. Expected Value of CHECKING with AA:")
    print("   - If Hero checks, the perfect Villain will make a minimum value bet with KK.")
    print("   - Hero (with AA) will raise, and Villain will fold.")
    print(f"   - Hero wins the pot. EV(Check) = ${ev_check_aa}")
    print("---------------------------------")

    # --- Step 2: Calculate EV of betting with AA ---
    # The EV of a value bet depends on the bet size S. A larger S is better if supported.
    # The derivative of EV(S) w.r.t S is positive, so we bet the max size possible.
    optimal_bet_size = stack

    # To make a bet credible, it must be balanced with bluffs. Villain must be indifferent
    # to calling the bet.
    # Villain's calling frequency 'c' is set to make Hero's bluffs break even.
    # EV(bluff) = c*(-S) + (1-c)*P = 0  => c*S = (1-c)*P => c*S = P - c*P => c(S+P) = P
    # c = P / (P + S)
    c = pot / (pot + optimal_bet_size)
    
    # EV for a value bet is EV(Bet) = c*(P+S) + (1-c)*P
    ev_bet_aa = c * (pot + optimal_bet_size) + (1 - c) * pot

    print("2. Expected Value of BETTING with AA:")
    print("   - The EV of value betting increases with bet size. So, the optimal bet size is all-in.")
    print(f"   - Optimal Bet Size (S) = ${optimal_bet_size}")
    print("   - To make this bet balanced, Villain must be made indifferent to calling.")
    print("   - This sets Villain's required calling frequency 'c'.")
    print(f"   - c = Pot / (Pot + S) = {pot} / ({pot} + {optimal_bet_size}) = {c:.4f}")
    print("   - The EV of the bet is then calculated:")
    # The final equation output as requested
    print(f"   - EV(Bet) = c * (Pot + S) + (1 - c) * Pot")
    print(f"   - EV(Bet) = {c:.4f} * ({pot} + {optimal_bet_size}) + (1 - {c:.4f}) * {pot} = ${ev_bet_aa:.2f}")
    print("---------------------------------")
    
    # --- Step 3: Compare EVs and determine strategy ---
    print("3. Strategy Derivation:")
    print(f"   - Comparing EVs: EV(Bet) = ${ev_bet_aa:.2f} vs. EV(Check) = ${ev_check_aa}")

    if ev_bet_aa > ev_check_aa:
        strategy_action = "BET"
        strategy_sizing = optimal_bet_size
        strategy_freq = 100
        print("   - Since betting has a higher EV, the optimal strategy for AA is to always bet.")
    elif ev_check_aa > ev_bet_aa:
        strategy_action = "CHECK"
        strategy_sizing = 0
        strategy_freq = 100
        print("   - Since checking has a higher EV, the optimal strategy for AA is to always check.")
    else:
        # This case implies indifference, but our calculation shows a clear preference.
        strategy_action = "BET" # Default to aggressive action
        strategy_sizing = optimal_bet_size
        strategy_freq = 100
        print("   - EVs are equal, can mix. But we default to the aggressive line.")

    # Rounding to nearest even number as per instructions
    rounded_sizing = int(round(strategy_sizing / 2) * 2)
    rounded_freq = int(round(strategy_freq / 2) * 2)
    
    print("---------------------------------")
    print("Final Strategy for AA:")
    final_strategy_str = ""
    if strategy_action == "BET":
        final_strategy_str = f"{strategy_action} {rounded_sizing} {rounded_freq}%"
    elif strategy_action == "CHECK":
        final_strategy_str = f"{strategy_action} {rounded_freq}%"

    print(f"   - Action: {strategy_action}, Sizing: {strategy_sizing}, Frequency: {strategy_freq}%")
    print(f"   - After rounding to nearest even number: Sizing={rounded_sizing}, Frequency={rounded_freq}%")
    print("\n<<<{}>>>".format(final_strategy_str))

solve_poker_strategy()
