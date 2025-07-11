import math

def solve_poker_scenario():
    """
    Calculates the GTO strategy for the given poker problem.
    """
    pot = 10
    stack = 1000
    
    # The optimal strategy is to bet the maximum possible amount (all-in).
    # This is found by maximizing Hero's EV function: EV = 10 - 50 / (Bet + 10)
    # The function is maximized when Bet is maximized.
    bet_size = stack
    
    # 1. Hero's strategy with AA (Value Hand)
    # With the nuts, Hero always wants to bet for value.
    hero_AA_bet_freq = 1.0 # 100%
    
    # 2. Hero's strategy with QQ (Bluff Hand)
    # Villain must be indifferent to calling. The bluff-to-value ratio must satisfy:
    # bluff_freq / value_freq = Bet / (Pot + Bet)
    # Since value_freq (p_A) is 1, bluff_freq (p_Q) = Bet / (Pot + Bet)
    hero_QQ_bet_freq_exact = bet_size / (pot + bet_size)
    hero_QQ_check_freq_exact = 1 - hero_QQ_bet_freq_exact
    
    # 3. Villain's response
    # Hero must be indifferent to bluffing. This determines Villain's calling frequency (MDF).
    # Villain_call_freq * (-Bet) + Villain_fold_freq * (Pot) = 0 (EV of bluffing)
    # c*(-B) + (1-c)*P = 0 => c = P / (P + B)
    villain_call_freq_exact = pot / (pot + bet_size)
    villain_fold_freq_exact = 1 - villain_call_freq_exact

    # Helper function for rounding to the nearest even number.
    # We find the two nearest even integers and see which one the number is closer to.
    def round_to_nearest_even(n):
        lower_even = math.floor(n / 2) * 2
        upper_even = math.ceil(n / 2) * 2
        if abs(n - lower_even) < abs(n - upper_even):
            return int(lower_even)
        elif abs(n - upper_even) < abs(n - lower_even):
            return int(upper_even)
        else:
            # In case of a tie (e.g., an odd integer), round half to the even number with an even half
            if (lower_even / 2) % 2 == 0:
                return int(lower_even)
            else:
                return int(upper_even)

    # A) Format Hero's action with AA
    hero_AA_bet_pct = round_to_nearest_even(hero_AA_bet_freq * 100)
    answer_A = f"BET ${int(bet_size)} {hero_AA_bet_pct}%"
    
    # B) Format Hero's action with QQ
    hero_QQ_bet_pct_raw = hero_QQ_bet_freq_exact * 100
    hero_QQ_bet_pct = round_to_nearest_even(hero_QQ_bet_pct_raw)
    hero_QQ_check_pct = 100 - hero_QQ_bet_pct
    
    # Sorted by amount of money put in (descending)
    answer_B_parts = []
    if hero_QQ_bet_pct > 0:
        answer_B_parts.append(f"BET ${int(bet_size)} {hero_QQ_bet_pct}%")
    if hero_QQ_check_pct > 0:
        answer_B_parts.append(f"CHECK {hero_QQ_check_pct}%")
    answer_B = ", ".join(answer_B_parts)

    # C) Format Villain's response
    villain_call_pct_raw = villain_call_freq_exact * 100
    villain_call_pct = round_to_nearest_even(villain_call_pct_raw)
    villain_fold_pct = 100 - villain_call_pct
    
    # Sorted by action (CALL then FOLD)
    answer_C_parts = []
    if villain_call_pct > 0:
        answer_C_parts.append(f"CALL {villain_call_pct}%")
    if villain_fold_pct > 0:
         answer_C_parts.append(f"FOLD {villain_fold_pct}%")
    answer_C = ", ".join(answer_C_parts)

    # Print the final formatted answer
    final_output = f"A) {answer_A} B) {answer_B} C) {answer_C}"
    print(final_output)

solve_poker_scenario()