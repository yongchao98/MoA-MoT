import math

def solve_poker_strategy():
    """
    Calculates the optimal GTO strategy for the given poker scenario.
    """
    pot = 10
    stack = 1000
    
    # We want to find the largest bet_size B that satisfies the constraints.
    # The EV-maximizing strategy is to use the largest possible bet size B.
    # A larger B corresponds to a larger bluffing frequency q, where q = B / (Pot + B).
    # We will iterate downwards from the highest possible even bluffing percentage (98%)
    # to find the first one that yields an even integer bet size.

    optimal_b = 0
    optimal_q_perc = 0

    # Iterate on q (bluff frequency with QQ) as an even percentage, from high to low.
    for q_perc in range(98, 0, -2):
        q = q_perc / 100.0
        
        # From the indifference principle q = B / (Pot + B), we derive B.
        # q * (Pot + B) = B  => q*Pot + q*B = B => q*Pot = B(1-q) => B = (q*Pot)/(1-q)
        b = (q * pot) / (1 - q)
        
        # Check if B is a whole number and an even number.
        if b.is_integer() and int(b) % 2 == 0:
            # Since we are iterating from the highest q, the first match
            # will be the one with the highest B, maximizing our EV.
            optimal_b = int(b)
            optimal_q_perc = q_perc
            break
            
    # Now we have the optimal bet size B and bluffing frequency q.
    # Hero's strategy (A, B) is determined.
    # A) With AA: Always bet for value. BET B 100%
    # B) With QQ: Mix between bluffing and checking.
    #    Bluff (BET B) with frequency q, check with frequency 1-q.
    
    hero_qq_check_perc = 100 - optimal_q_perc
    
    # C) Villain's Response
    # To make Hero indifferent to bluffing with QQ, Villain must mix calling/folding.
    # EV(bluff) = c*(-B) + (1-c)*Pot = 0 (where c is call frequency)
    # c*(-B) + Pot - c*Pot = 0 => Pot = c*(B + Pot) => c = Pot / (B + Pot)
    
    villain_call_freq = pot / (optimal_b + pot)
    villain_call_perc = round(villain_call_freq * 100)
    villain_fold_perc = 100 - villain_call_perc
    
    # Format the final answer
    # A) you hold AA
    # B) you hold QQ
    # C) villains response to your bets
    
    # Using RAISE, BET, CALL, CHECK and FOLD for the actions.
    # Sort mixed strategies by money put in (CHECK then BET).
    
    print(f"A) BET ${optimal_b} 100%")
    print(f"B) CHECK {hero_qq_check_perc}%, BET ${optimal_b} {optimal_q_perc}%")
    print(f"C) CALL {villain_call_perc}%, FOLD {villain_fold_perc}%")
    
    # Final check of numbers in the equation:
    # B = (q*Pot)/(1-q) => 240 = (0.96 * 10) / (1 - 0.96) = 9.6 / 0.04 = 240. Correct.
    # c = Pot / (B + Pot) => c = 10 / (240 + 10) = 10 / 250 = 0.04 = 4%. Correct.

solve_poker_strategy()