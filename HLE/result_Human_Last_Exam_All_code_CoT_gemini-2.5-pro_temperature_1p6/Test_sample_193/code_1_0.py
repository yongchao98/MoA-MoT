import math

def solve_poker_strategy():
    """
    Analyzes the poker scenario to find the optimal strategy for holding AA.
    """
    # Game parameters
    POT_SIZE = 10
    STACK_SIZE = 1000
    
    # Range analysis
    # Hero has 1 value hand type (AA) and 10 bluff hand types (QQ-33)
    # vs Villain's KK on a 222 board.
    
    # --- Step 1: Analyze EV of Checking with AA ---
    # If we check, we go to showdown and win the pot.
    ev_check_aa = POT_SIZE
    
    print("--- Analysis for Hero holding AA ---")
    print("Option 1: CHECK")
    print(f"If Hero checks with AA, they win the pot at showdown.")
    print(f"The equation for the Expected Value of checking is:")
    print(f"EV(CHECK) = Pot = {ev_check_aa}\n")

    # --- Step 2: Analyze EV of Betting with AA ---
    # In a GTO strategy, a larger bet size for a value hand leads to higher EV.
    # Therefore, the optimal bet size is the largest possible: all-in.
    bet_size = STACK_SIZE

    # To be unexploitable, this bet must be balanced with bluffs. Villain's optimal
    # counter-strategy is to call with a frequency that makes our bluffs break-even.
    # Villain's Call Frequency (c) is such that EV(Hero Bluff) = 0.
    # c * (-Bet) + (1-c) * (Pot) = 0  => c * Bet = (1-c) * Pot => c * (Bet + Pot) = Pot
    # c = Pot / (Pot + Bet)
    villain_call_freq = POT_SIZE / (POT_SIZE + bet_size)
    villain_fold_freq = 1 - villain_call_freq

    # Now, calculate the EV of betting all-in with AA given Villain's optimal calling.
    ev_win_on_call = POT_SIZE + bet_size
    ev_win_on_fold = POT_SIZE
    
    ev_bet_aa = (villain_call_freq * ev_win_on_call) + (villain_fold_freq * ev_win_on_fold)

    print("Option 2: BET")
    print(f"To maximize value, Hero should consider betting the largest possible amount, which is all-in for {bet_size}.")
    print("The equation for the Expected Value of betting is based on Villain's optimal calling frequency:")
    print(f"EV(BET) = (Villain Calls % * Pot Won on Call) + (Villain Folds % * Pot Won on Fold)")
    # Output each number in the final equation
    print(f"EV(BET {bet_size}) = ({villain_call_freq:.4f} * {ev_win_on_call}) + ({villain_fold_freq:.4f} * {ev_win_on_fold})")
    print(f"EV(BET {bet_size}) = {ev_bet_aa:.2f}\n")
    
    # --- Step 3: Compare EVs and determine strategy ---
    print("--- Conclusion ---")
    print(f"Comparing the Expected Values: EV(BET) of ${ev_bet_aa:.2f} is much greater than EV(CHECK) of ${ev_check_aa}.")
    print("Therefore, the optimal strategy with AA is always to bet all-in.\n")

    action = "BET"
    sizing = bet_size
    frequency_percent = 100

    # Round sizing and frequency to the nearest even number per instructions.
    final_sizing = int(round(sizing / 2.0) * 2)
    final_freq = int(round(frequency_percent / 2.0) * 2)

    print("--- Final Formatted Answer ---")
    print(f"{action} {final_sizing} {final_freq}%")


solve_poker_strategy()