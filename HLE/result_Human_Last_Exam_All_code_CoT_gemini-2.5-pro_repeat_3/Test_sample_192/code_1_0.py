import math

def solve_poker_strategy():
    """
    Calculates and explains the optimal poker strategy for the given river scenario.
    """
    pot = 10
    stack = 1000

    # Step 1: Determine the optimal bet size.
    # To maximize value from our winning hands (AA), we should bet as large as possible.
    bet_size = stack
    print(f"Step 1: Determine the optimal bet size.")
    print(f"To maximize value from our winning hands (AA), we bet the entire stack.")
    print(f"Optimal Bet Size = ${bet_size}\n")

    # Step 2: Calculate the GTO bluffing frequency for the losing hand (QQ).
    # The formula to make Villain indifferent is: Frequency = Bet / (Pot + Bet)
    # Since we always bet our value hand (AA), this dictates the bluffing frequency for QQ.
    print(f"Step 2: Calculate the GTO bluffing frequency for QQ.")
    print(f"The formula for the bluffing frequency is Bet / (Pot + Bet).")
    bluff_freq_raw = bet_size / (pot + bet_size)
    print(f"Frequency = {bet_size} / ({pot} + {bet_size}) = {bluff_freq_raw:.4f}\n")

    bluff_pct_raw = bluff_freq_raw * 100
    print(f"Step 3: Apply the rounding constraint ('nearest even number').")
    print(f"The ideal bluffing percentage is {bluff_pct_raw:.2f}%.")
    
    # The nearest even integers to 99.01 are 98 and 100.
    # We must analyze the outcome of both choices against a 'perfect nemesis'.
    pct_up = 100
    pct_down = 98

    print(f"The nearest even integers are {pct_down}% and {pct_up}%. We must choose the one with the higher Expected Value (EV).\n")

    print(f"Analysis for rounding UP to {pct_up}%:")
    print(f"If we bluff {pct_up}% of the time, our bluff-to-value ratio is {pct_up/100:.2f}, which is higher than the GTO ratio of {bluff_freq_raw:.2f} (over-bluffing).")
    print(f"A perfect nemesis will exploit this by always CALLING.")
    # EV = 0.5 * (EV with AA) + 0.5 * (EV with QQ) = 0.5 * (Pot + Bet) + 0.5 * (-Bet) = 0.5 * Pot
    ev_up = 0.5 * pot
    print(f"Our EV = 0.5 * (${pot+bet_size}) + 0.5 * (-${bet_size}) = ${ev_up:.2f}\n")

    print(f"Analysis for rounding DOWN to {pct_down}%:")
    print(f"If we bluff {pct_down}% of the time, our ratio is {pct_down/100:.2f}, which is lower than the GTO ratio of {bluff_freq_raw:.2f} (under-bluffing).")
    print(f"A perfect nemesis will exploit this by always FOLDING to our bet.")
    # EV_AA = Pot (villain folds)
    # EV_QQ = 0.98 * Pot (bet, villain folds) + 0.02 * 0 (check, lose pot)
    ev_aa_down = pot
    ev_qq_down = (pct_down / 100) * pot
    ev_down = 0.5 * ev_aa_down + 0.5 * ev_qq_down
    print(f"Our EV = 0.5 * (${ev_aa_down}) + 0.5 * (({pct_down}/100) * ${pot}) = ${ev_down:.2f}\n")

    print(f"Step 4: Conclude the optimal strategy.")
    print(f"Comparing the EVs, ${ev_down:.2f} (from rounding down) is much better than ${ev_up:.2f} (from rounding up).")
    print(f"Therefore, the optimal strategy is to bluff with QQ at {pct_down}% frequency, which forces the villain to fold.\n")

    # Final strategy based on the higher EV choice
    bet_aa_pct = 100
    bet_qq_pct = pct_down
    check_qq_pct = 100 - bet_qq_pct
    villain_response = "FOLD 100%"

    # Format the final answer string
    ans_A = f"BET ${bet_size} {bet_aa_pct}%"
    ans_B = f"CHECK {check_qq_pct}%, BET ${bet_size} {bet_qq_pct}%"
    ans_C = villain_response

    print("----------------------------------------")
    print("Final Answer:")
    print(f"A) {ans_A}")
    print(f"B) {ans_B}")
    print(f"C) {ans_C}")
    print("----------------------------------------")

solve_poker_strategy()