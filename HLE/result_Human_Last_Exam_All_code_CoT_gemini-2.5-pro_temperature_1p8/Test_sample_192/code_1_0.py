def solve_poker_problem():
    """
    Calculates the optimal GTO strategy for the given poker scenario,
    adhering to the even number rounding constraint.
    """
    pot_size = 10

    # To maximize EV, we need to maximize bet size B.
    # The constraint is that all frequencies must be even integers.
    # Our bluffing frequency with QQ is p = B / (pot + B).
    # To maximize B, we must maximize p. The highest possible even
    # integer percentage for p below 100% is 98%.
    p_bluff_frequency = 0.98 # 98%

    # From p = B / (pot + B), we solve for B:
    # p * (pot + B) = B
    # p*pot + p*B = B
    # p*pot = B - p*B
    # p*pot = B * (1 - p)
    # B = (p * pot) / (1 - p)
    bet_size = (p_bluff_frequency * pot_size) / (1 - p_bluff_frequency)

    # Villain's calling frequency 'q' is what makes us indifferent to bluffing.
    # q = pot / (pot + B)
    q_call_frequency = pot_size / (pot_size + bet_size)

    # Sanity check our values
    # bet_size should be 490.
    # q_call_frequency should be 10 / 500 = 0.02
    
    print(f"To satisfy the even-number constraint while maximizing EV, we find the highest even bluffing percentage possible, which is {int(p_bluff_frequency*100)}%.")
    print("This leads to the following equation for the bet size (B):")
    print(f"B = ({p_bluff_frequency} * {pot_size}) / (1 - {p_bluff_frequency}) = {int(bet_size)}")
    print("\nThe villain's corresponding optimal calling frequency (q) is:")
    print(f"q = {pot_size} / ({pot_size} + {int(bet_size)}) = {q_call_frequency}")
    
    # Format the final answer string
    
    # A) With AA, we always bet for value. 100% is an even number.
    action_A = f"BET ${int(bet_size)} 100%"

    # B) With QQ, we bluff with frequency p and check otherwise.
    # The actions are sorted by money put in (Bet > Check).
    action_B = f"BET ${int(bet_size)} {int(p_bluff_frequency*100)}%, CHECK {int((1-p_bluff_frequency)*100)}%"

    # C) Villain's response. They call with frequency q and fold otherwise.
    # Sorted by money put in (Call > Fold).
    action_C = f"CALL {int(q_call_frequency*100)}%, FOLD {int((1-q_call_frequency)*100)}%"
    
    final_answer = f"A) {action_A} B) {action_B} C) {action_C}"

    print("\nFinal Answer:")
    print(final_answer)

solve_poker_problem()
<<<A) BET $490 100% B) BET $490 98%, CHECK 2% C) CALL 2%, FOLD 98%>>>