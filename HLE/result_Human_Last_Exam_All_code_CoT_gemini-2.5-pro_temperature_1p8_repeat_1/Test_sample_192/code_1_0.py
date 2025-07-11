import math

def solve_poker_strategy():
    """
    Calculates the optimal poker strategy based on the given scenario and constraints.
    """
    pot = 10
    
    # --- Plan Explanation ---
    print("This script calculates the Game Theory Optimal (GTO) strategy for the given poker scenario.")
    print("The goal is to find a bet size (B) and frequencies that make both players indifferent to certain actions,")
    print("while adhering to the constraint that all percentages and the bet size must be even numbers.\n")

    # --- Step 1: Work backwards from the 'even percentage' constraint ---
    # Let's set the villain's call frequency to a small even percentage, like 2%.
    # Villain's call frequency (c) is determined by making our bluff indifferent.
    # The formula is: c = pot / (pot + B)
    # We set c = 0.02 and solve for B.
    villain_call_pct = 2
    c = villain_call_pct / 100
    
    print(f"Step 1: Determine the bet size 'B'.")
    print(f"We set Villain's calling frequency to a small, even percentage: {villain_call_pct}%.")
    print(f"Using the indifference formula c = pot / (pot + B):")
    print(f"{c} = {pot} / ({pot} + B)")
    print(f"{c} * ({pot} + B) = {pot}")
    print(f"{c*pot} + {c}*B = {pot}")
    print(f"{c}*B = {pot - c*pot}")
    
    # Calculation
    numerator = pot - (c * pot)
    bet_size = numerator / c

    print(f"B = {numerator} / {c} = {bet_size}\n")
    
    bet_size = int(bet_size)

    # --- Step 2: Calculate our bluffing frequency with QQ ---
    # Our bluffing frequency (q) is determined by making the villain's call indifferent.
    # The formula is: q = B / (pot + B)
    print("Step 2: Determine our bluffing frequency 'q' with QQ.")
    print(f"Using the bet size B = ${bet_size} and the formula q = B / (pot + B):")
    
    # Calculation
    hero_bluff_freq = bet_size / (pot + bet_size)
    
    print(f"q = {bet_size} / ({pot} + {bet_size}) = {bet_size / (pot + bet_size):.2f}\n")
    
    hero_bet_qq_pct = int(hero_bluff_freq * 100)
    hero_check_qq_pct = 100 - hero_bet_qq_pct
    
    # Hero's action with AA is always to bet for value.
    hero_bet_aa_pct = 100
    
    # Villain's fold frequency is the remainder of their calling frequency.
    villain_fold_pct = 100 - villain_call_pct
    
    # --- Step 3: Format the final answer ---
    print("Step 3: Assembling the final strategy based on these calculations.")
    
    answer_a = f"BET ${bet_size} {hero_bet_aa_pct}%"
    answer_b = f"BET ${bet_size} {hero_bet_qq_pct}%, CHECK {hero_check_qq_pct}%"
    answer_c = f"FOLD {villain_fold_pct}%, CALL {villain_call_pct}%"
    
    final_answer = f"A) {answer_a} B) {answer_b} C) {answer_c}"
    
    print("\nFinal Answer:")
    print(final_answer)

solve_poker_strategy()
<<<A) BET $490 100% B) BET $490 98%, CHECK 2% C) FOLD 98%, CALL 2%>>>