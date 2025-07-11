def analyze_shove_options():
    """
    Analyzes poker hand options for an all-in shove near the money bubble.
    """
    # --- Situation ---
    position = "UTG+1"
    stack_bb = 16
    situation = "Near the money bubble (high ICM pressure)"

    # --- Analysis ---
    # On the bubble, opponents will call an early position shove with a very tight range.
    # A reasonable assumption for a calling range is TT+, AQs+, AKo.
    # We need to shove a hand that has good equity WHEN CALLED by this range.
    calling_range = "TT+, AQs+, AKo"

    # Hand options and their approximate equity against the calling range.
    # Equity is the percentage chance to win the pot if called.
    hand_options = {
        "QJs": 33,
        "99": 36,
        "AJo": 28,
        "AKo": 47
    }

    print(f"Analyzing shove from {position} with {stack_bb}bb stack.")
    print(f"Situation: {situation}")
    print(f"Assumed opponent calling range: {calling_range}\n")
    print("--- Evaluating Hand Options ---")

    # A. QJs
    hand_a = "QJs"
    equity_a = hand_options[hand_a]
    print(f"A. {hand_a}: Equity vs calling range is ~{equity_a}%.")
    print("   This is too low. It's a speculative hand that performs poorly when called. FOLD.\n")

    # C. 99
    hand_c = "99"
    equity_c = hand_options[hand_c]
    print(f"C. {hand_c}: Equity vs calling range is ~{equity_c}%.")
    print("   While a strong pair, it is behind all pairs in the calling range (TT, JJ, QQ, KK, AA).")
    print("   This shove relies heavily on opponents folding. RISKY.\n")

    # D. AJo
    hand_d = "AJo"
    equity_d = hand_options[hand_d]
    print(f"D. {hand_d}: Equity vs calling range is ~{equity_d}%.")
    print("   This hand is dominated by AQ and AK in the calling range. VERY BAD SHOVE.\n")

    # E. AKo
    hand_e = "AKo"
    equity_e = hand_options[hand_e]
    print(f"E. {hand_e}: Equity vs calling range is ~{equity_e}%.")
    print("   This is nearly a coin flip against the entire calling range.")
    print("   It also blocks Aces and Kings, making it less likely opponents have a premium hand.")
    print("   This is the most robust and profitable shove.\n")

    print("--- Conclusion ---")
    print("AKo is the strongest option. It has the best equity against a tight calling range and good blocker value, making it the correct hand to jam all-in.")

analyze_shove_options()
<<<E>>>