def solve_poker_strategy():
    """
    This function analyzes the poker scenario to determine the optimal strategy for AA.
    It prints the logical steps and the final conclusion.
    """

    # --- Game State ---
    pot = 10
    villain_bet_if_checked = 0.02 # A minimal, non-zero bet (epsilon)

    # --- Analysis for AA ---

    # 1. Calculate the EV of Checking with AA
    ev_check_aa = pot + villain_bet_if_checked
    print("Step 1: Calculate EV of checking with AA.")
    print(f"If we check, Villain's optimal response is to make a tiny bet ({villain_bet_if_checked}). We call and win the pot plus this bet.")
    print(f"EV_check(AA) = {pot} + {villain_bet_if_checked} = {ev_check_aa}\n")

    # 2. Calculate the EV of Betting with AA
    # To maximize EV from betting, we should use the largest possible bet size
    # that we can balance with our available bluffs.
    # We have 10 bluff combos (QQ-33) to balance 1 value combo (AA).
    # Number of bluffs required (Rb) = Bet / Pot
    # Maximum Rb is 10, so max Bet = 10 * Pot = 100.
    optimal_bet_size = 100
    
    # The formula for EV of betting a value hand in this spot is:
    # EV_bet(AA) = Pot + (Pot * Bet) / (Pot + Bet)
    ev_bet_aa = pot + (pot * optimal_bet_size) / (pot + optimal_bet_size)

    print("Step 2: Calculate EV of betting with AA.")
    print("Betting is better with larger sizes. The max supportable bet size is 100.")
    print(f"The calculation for EV_bet(AA) with bet size B={optimal_bet_size} is:")
    print(f"EV_bet(AA) = {pot} + ({pot} * {optimal_bet_size}) / ({pot} + {optimal_bet_size})")
    print(f"EV_bet(AA) = {ev_bet_aa:.2f}\n")

    # 3. Compare EVs and Conclude
    print("Step 3: Compare EVs and determine the optimal strategy.")
    print(f"Comparing EV_bet({ev_bet_aa:.2f}) with EV_check({ev_check_aa:.2f}).")
    if ev_bet_aa > ev_check_aa:
        print("Betting is strictly more profitable than checking.")
        action = "BET"
        sizing = optimal_bet_size
        percentage = 100
    else:
        # This case is not reached with this logic, but included for completeness.
        print("Checking is the optimal action.")
        action = "CHECK"
        sizing = 0
        percentage = 100

    # Rounding to nearest even number as per instruction.
    # 100 is already an even number for both sizing and percentage.
    final_sizing = int(round(sizing / 2.0) * 2)
    final_percentage = int(round(percentage / 2.0) * 2)

    print("\nFinal optimal strategy for AA is:")
    print(f"{action} {final_sizing} {final_percentage}%")
    
    # Store the final answer in the required format for the very end.
    global final_answer
    final_answer = f"{action} {final_sizing} {final_percentage}%"


# Run the analysis
solve_poker_strategy()
print(f"\n<<< {final_answer} >>>")