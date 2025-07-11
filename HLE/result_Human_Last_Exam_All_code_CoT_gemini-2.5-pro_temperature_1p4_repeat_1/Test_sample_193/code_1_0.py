import math

def solve_poker_strategy():
    """
    Calculates the optimal poker strategy for the given scenario.
    """
    
    # --- Step 1 & 2: Define the Game State and Problem ---
    pot_size = 10
    stack_size = 1000
    hero_value_combos = 1  # AA is the only hand that beats KK
    hero_bluff_combos = 10 # QQ, JJ, TT, 99, 88, 77, 66, 55, 44, 33

    print("--- Step-by-step GTO Analysis ---")
    print(f"Pot Size (P): {pot_size}")
    print(f"Effective Stack Size: {stack_size}")
    print(f"Hero's Value Combos (V): {hero_value_combos} (AA)")
    print(f"Hero's Bluff Candidates: {hero_bluff_combos} (QQ-33)")

    # --- Step 3: Villain's Indifference Equation ---
    print("\nStep 3: To make Villain indifferent to calling our bet (B), we must balance our betting range.")
    print("The required ratio of Bluffs (Bl) to Value (V) is:")
    print("Bl / V = B / (P + B)")
    print(f"Since V = {hero_value_combos} and P = {pot_size}, the number of bluff combos we must bet is:")
    print(f"Bl = {hero_value_combos} * B / ({pot_size} + B) => Bl = B / ({pot_size} + B)")

    # --- Step 4: Maximizing Hero's EV ---
    print("\nStep 4: We choose a bet size B to maximize our range's Expected Value (EV).")
    print("Assuming Villain calls when indifferent, the EV of our range is:")
    print("EV(B) = (EV from Value Hands) + (EV from Bluffs)")
    print(f"EV(B) = V*(P + B) + Bl*(-B)")
    print(f"EV(B) = {hero_value_combos}*({pot_size} + B) - (B / ({pot_size} + B)) * B")
    print(f"EV(B) = {pot_size} + B - B^2 / ({pot_size} + B)")
    
    print("\nTo maximize EV(B), we take its derivative with respect to B and analyze it:")
    # The derivative is d(EV)/dB = P^2 / (P + B)^2
    derivative_numerator = pot_size**2
    print(f"d(EV)/dB = {derivative_numerator} / ({pot_size} + B)^2")

    print("\nSince the derivative is always positive for any B > 0, the function EV(B) is always increasing.")
    print("Therefore, to maximize our EV, we must choose the largest possible bet size.")
    
    # --- Step 5: Determine the Final Strategy ---
    optimal_bet_size = stack_size # Bet all-in
    
    print(f"\nThe largest possible bet is the entire stack: {optimal_bet_size}.")
    print("This means the optimal strategy with our value hand (AA) is to bet the maximum amount to extract the most value.")
    
    # Calculate strategy for AA
    action = "BET"
    sizing = optimal_bet_size
    frequency = 100
    
    # Round to nearest even number as per instructions
    final_sizing = int(2 * round(sizing/2))
    final_frequency = int(2 * round(frequency/2))
    
    print("\n--- Optimal Strategy for AA ---")
    print(f"The analysis shows a pure strategy for AA. We do not need to mix actions.")
    print("Action: BET")
    print(f"Sizing: {final_sizing}")
    print(f"Frequency: {final_frequency}%")
    
    # Final answer in the required format
    final_answer = f"{action} {final_sizing} {final_frequency}%"
    print(f"\n<<<"+final_answer+">>>")

solve_poker_strategy()