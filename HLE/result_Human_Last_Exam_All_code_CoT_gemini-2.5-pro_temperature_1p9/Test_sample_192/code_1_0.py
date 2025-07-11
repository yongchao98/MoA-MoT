import math

def solve_poker_problem():
    """
    Calculates the optimal poker strategy based on the given constraints.

    The logic finds a bet size `B` that results in GTO frequencies that are
    naturally even percentages, satisfying the problem's rounding constraint
    and maintaining a perfect equilibrium.
    """
    pot = 10.0
    # We need to find a bet_size B such that all frequencies are even integers.
    # Villain's calling frequency F_call = pot / (pot + B)
    # We test for the smallest non-zero even percentage, 2% (0.02), for F_call.
    # 0.02 = 10 / (10 + B) => 0.2 + 0.02B = 10 => 0.02B = 9.8 => B = 490
    
    bet_size = 490

    # Case A: You hold AA (Value Bet)
    # You always bet your nut hands for value.
    action_A = f"BET ${int(bet_size)} 100%"

    # Case B: You hold QQ (Bluff or Check)
    # Your bluffing frequency f_bluff = bet_size / (pot + bet_size)
    # Your checking frequency f_check = 1 - f_bluff
    f_bluff = bet_size / (pot + bet_size)
    bluff_pct = int(f_bluff * 100)
    check_pct = 100 - bluff_pct
    
    # Sort by amount put in (CHECK is $0, BET is > $0)
    action_B = f"CHECK {check_pct}%, BET ${int(bet_size)} {bluff_pct}%"

    # Case C: Villain's response
    # Villain's calling frequency F_call = pot / (pot + bet_size)
    # Villain's folding frequency F_fold = 1 - F_call
    f_call = pot / (pot + bet_size)
    call_pct = int(f_call * 100)
    fold_pct = 100 - call_pct
    action_C = f"CALL {call_pct}%, FOLD {fold_pct}%"

    final_answer = f"A) {action_A} B) {action_B} C) {action_C}"
    
    print(final_answer)

solve_poker_problem()