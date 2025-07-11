import math

def solve_poker_strategy():
    """
    Calculates the optimal GTO strategy for the given poker scenario.
    """

    # Step 1: Define problem parameters
    P = 10.0          # Pot size
    stack = 1000.0    # Effective stack size
    V = 1.0           # Number of value combos (AA)
    L_total = 10.0    # Number of total bluff combos (QQ-33)

    print("--- GTO Poker Strategy Calculation ---")
    print(f"Pot: ${int(P)}, Stack: ${int(stack)}, Hero Value Combos: {int(V)}, Hero Bluff Combos: {int(L_total)}")
    print("-" * 38)

    # Step 2: Explain the theory
    print("\n[Analysis]")
    print("Hero has a polarized range (nuts or bluffs) vs. Villain's bluff-catcher.")
    print("The optimal strategy is to bet a size 'B' with all value hands and a number of bluffs 'L' to make Villain indifferent to calling.")
    print("The required bluff-to-value ratio is: L/V = B / (P + B)")

    # Step 3: Explain EV Maximization
    print("\nThe EV of this strategy increases as the bet size 'B' increases because the derivative dEV/dB = P^2 / (P + B)^2 is always positive.")
    print("Therefore, the optimal strategy is to bet the largest amount possible: the entire stack.")

    # Step 4: Determine optimal bet size and check constraints
    optimal_bet_size = stack
    print(f"\nOptimal Bet Size (B) = ${int(optimal_bet_size)}")

    print("\n[Constraint Check]")
    print("We must verify we have enough bluffs to balance this bet size.")
    # The final equation for the number of bluffs needed (L)
    L_needed = V * optimal_bet_size / (P + optimal_bet_size)
    
    print(f"Equation for bluffs needed: L = V * B / (P + B)")
    print(f"Calculation: L = {int(V)} * {int(optimal_bet_size)} / ({int(P)} + {int(optimal_bet_size)}) = {L_needed:.2f}")

    if L_needed <= L_total:
        print(f"Bluffs needed ({L_needed:.2f}) is less than or equal to available bluffs ({int(L_total)}). The strategy is valid.")
    else:
        print("Error: Not enough bluffs to support this bet size.")
        return

    # Step 5: Format the final answer
    action = "BET"
    # Round sizing and percentage to the nearest even number
    final_sizing = int(round(optimal_bet_size / 2) * 2)
    final_frequency = int(round(100 / 2) * 2)

    print("\n--- Final Optimal Strategy for AA ---")
    print("Action: BET")
    print(f"Sizing: {final_sizing}")
    print(f"Frequency: {final_frequency}%")

    final_answer = f"{action} {final_sizing} {final_frequency}%"
    print(f"\nFormatted Answer: {final_answer}")

solve_poker_strategy()