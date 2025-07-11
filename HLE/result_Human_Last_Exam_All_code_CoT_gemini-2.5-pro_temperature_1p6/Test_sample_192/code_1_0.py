import math

def round_to_nearest_even(n):
    """Rounds a number to the nearest even integer."""
    return int(round(n / 2.0) * 2.0)

def solve_poker_strategy():
    """
    Calculates the GTO strategy for the given poker scenario and prints the result.
    """
    pot = 10
    stack = 1000

    # 1. Optimal bet size is the effective stack to maximize EV.
    bet_size = stack

    # 2. Calculate Villain's calling frequency.
    # Villain must call with a frequency that makes Hero indifferent to bluffing.
    # EV(bluff) = EV(check) => p_call * (-bet_size) + (1 - p_call) * pot = 0
    # This solves to p_call = pot / (pot + bet_size)
    p_call = pot / (pot + bet_size)
    p_fold = 1 - p_call
    
    # Convert to percentages for formatting.
    p_call_pct = p_call * 100
    p_fold_pct = p_fold * 100
    
    # 3. Calculate Hero's bluffing frequency.
    # Hero must bluff with a frequency that makes Villain indifferent to calling.
    # Villain is indifferent when P(Hero is bluffing | Hero bets) = Villain's Pot Odds
    # P(bluff|bet) = bet_size / (pot + 2 * bet_size)
    # Since Hero's range is 50/50 and Hero always bets value (AA), we have:
    # p_bet_QQ / (1 + p_bet_QQ) = bet_size / (pot + 2 * bet_size)
    bluff_ratio_in_betting_range = bet_size / (pot + 2 * bet_size)
    # Solving for p_bet_QQ: p = r / (1 - r)
    p_bet_QQ = bluff_ratio_in_betting_range / (1 - bluff_ratio_in_betting_range)
    p_check_QQ = 1 - p_bet_QQ

    # Convert to percentages for formatting.
    p_bet_QQ_pct = p_bet_QQ * 100
    p_check_QQ_pct = p_check_QQ * 100

    # 4. Apply rounding rules.
    # To ensure percentages sum to 100, we round the larger value and derive the smaller.
    if p_bet_QQ_pct > p_check_QQ_pct:
        rounded_p_bet_QQ_pct = round_to_nearest_even(p_bet_QQ_pct)
        rounded_p_check_QQ_pct = 100 - rounded_p_bet_QQ_pct
    else:
        rounded_p_check_QQ_pct = round_to_nearest_even(p_check_QQ_pct)
        rounded_p_bet_QQ_pct = 100 - rounded_p_check_QQ_pct
        
    if p_call_pct > p_fold_pct:
        rounded_p_call_pct = round_to_nearest_even(p_call_pct)
        rounded_p_fold_pct = 100 - rounded_p_call_pct
    else:
        rounded_p_fold_pct = round_to_nearest_even(p_fold_pct)
        rounded_p_call_pct = 100 - rounded_p_fold_pct

    # 5. Format and print the final output.

    # A) Strategy with AA (always bet for value).
    ans_A = f"BET ${int(bet_size)} 100%"

    # B) Strategy with QQ (mixed). Sorted by money committed: BET ($1000) then CHECK ($0).
    ans_B = f"BET ${int(bet_size)} {rounded_p_bet_QQ_pct}%, CHECK {rounded_p_check_QQ_pct}%"

    # C) Villain's response. Sorted by money committed: CALL ($1000) then FOLD ($0).
    # Since call % is 0, we can show Fold first for clarity, but strictly following money sort:
    if rounded_p_call_pct > 0:
        ans_C = f"CALL {rounded_p_call_pct}%, FOLD {rounded_p_fold_pct}%"
    else:
        ans_C = f"FOLD {rounded_p_fold_pct}%, CALL {rounded_p_call_pct}%"


    final_output = f"""A) {ans_A}
B) {ans_B}
C) {ans_C}
<<<A) {ans_A} B) {ans_B} C) {ans_C}>>>"""
    
    print(final_output)

solve_poker_strategy()