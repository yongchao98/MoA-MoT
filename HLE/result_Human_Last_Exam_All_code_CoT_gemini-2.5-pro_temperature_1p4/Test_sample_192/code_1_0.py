import math

def solve_poker_strategy():
    """
    Solves the GTO poker problem based on the provided scenario and constraints.
    
    The solution is derived by finding an equilibrium that satisfies two key 
    indifference equations while adhering to the constraint that all action
    frequencies must be even percentages.
    
    1. Villain's calling frequency (f_call) makes Hero indifferent to bluffing with QQ.
       f_call = Pot / (Pot + B)
    
    2. Hero's bluffing frequency with QQ (p_QQ) makes Villain indifferent to calling Hero's bet.
       p_QQ = B / (Pot + B)

    By working backwards from the even-percentage constraint, we find the optimal
    bet size (B) and the corresponding frequencies for both players.
    """
    pot = 10.0
    optimal_B = 0
    p_bluff_pct = 0
    f_call_pct = 0

    # We iterate through possible even percentages for the villain's calling frequency
    # to find one that results in a valid, constrained equilibrium.
    for call_pct_int in range(2, 100, 2):
        f_call = call_pct_int / 100.0
        
        # From f_call = Pot / (Pot + B), we solve for B
        # B = (Pot / f_call) - Pot
        # Note: B must be > 0, so f_call cannot be >= 1.
        if f_call > 0:
            B = (pot / f_call) - pot
        else:
            continue

        # Check if B is an even integer (satisfying sizing constraint).
        if B > 0 and round(B) == B and round(B) % 2 == 0:
            # Check if this B gives an even bluffing percentage for Hero.
            # p_QQ = B / (Pot + B)
            p_bluff = B / (pot + B)
            p_bluff_pct_val = p_bluff * 100
            
            p_check_pct_val = (1 - p_bluff) * 100

            # Check if both bluffing and checking percentages are even integers.
            if round(p_bluff_pct_val) == p_bluff_pct_val and p_bluff_pct_val % 2 == 0:
                 if round(p_check_pct_val) == p_check_pct_val and p_check_pct_val % 2 == 0:
                    optimal_B = int(B)
                    p_bluff_pct = int(p_bluff_pct_val)
                    f_call_pct = call_pct_int
                    break # Solution found

    # Assemble the final answer string according to the specified format.
    # The actions are sorted by the amount of money put in, from largest to smallest.

    # A) With AA (the nuts), we always bet for value to balance our bluffs.
    ans_A = f"BET ${optimal_B} 100%"

    # B) With QQ, we bluff at the calculated frequency.
    p_check_pct = 100 - p_bluff_pct
    ans_B = f"BET ${optimal_B} {p_bluff_pct}%, CHECK {p_check_pct}%"

    # C) Villain's response is to call at the calculated frequency.
    f_fold_pct = 100 - f_call_pct
    ans_C = f"CALL {f_call_pct}%, FOLD {f_fold_pct}%"

    final_output = f"A) {ans_A} B) {ans_B} C) {ans_C}"
    print(final_output)

solve_poker_strategy()