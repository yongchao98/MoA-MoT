import sys

def solve_utility():
    """
    Solves the Alice and Bob game theory problem by applying
    Iterated Elimination of Weakly Dominated Strategies (IEWDS).
    """
    # Define the actions and the payoff matrix
    # Payoffs are structured as (Alice's Payoff, Bob's Payoff)
    actions = ["Rest", "Bike", "Run"]
    payoffs = [
        # Bob's actions: Rest, Bike, Run
        [(0, 0), (2, 0), (4, 0)],  # Alice Rests
        [(0, 2), (-2, -2), (2, 0)],  # Alice Bikes
        [(0, 4), (0, 2), (-3, -3)],  # Alice Runs
    ]

    print("--- Initial 3x3 Payoff Matrix (Alice, Bob) ---")
    print("           Bob: Rest   Bike     Run")
    for i, row in enumerate(payoffs):
        print(f"Alice {actions[i]:<4}: {str(row[0]):<8} {str(row[1]):<8} {str(row[2]):<8}")
    print("-" * 48)

    # Convert to a more usable dictionary format for elimination
    alice_available = {i: [payoffs[i][j][0] for j in range(3)] for i in range(3)}
    bob_available = {j: [payoffs[i][j][1] for i in range(3)] for j in range(3)}
    
    # --- Round 1 of Elimination ---
    print("Step 1: Analyzing the full 3x3 game for weakly dominated strategies.\n")
    
    # For Alice, Rest (0, 2, 4) weakly dominates Run (0, 0, -3)
    # because 0==0, 2>0, and 4>-3.
    # Alice, being superrational, will never choose 'Run'.
    print("For Alice: Comparing 'Rest' ([0, 2, 4]) vs 'Run' ([0, 0, -3])")
    print("'Rest' payoffs are equal or better than 'Run' payoffs in all cases, and strictly better in some.")
    print("=> Alice eliminates 'Run' as a choice.\n")
    alice_available.pop(2) # Remove Run

    # For Bob (symmetric game), Rest (0, 2, 4) weakly dominates Run (0, 0, -3).
    # Bob, also superrational, will never choose 'Run'.
    print("For Bob: The game is symmetric. Comparing 'Rest' ([0, 2, 4]) vs 'Run' ([0, 0, -3])")
    print("'Rest' payoffs are equal or better than 'Run' payoffs in all cases, and strictly better in some.")
    print("=> Bob eliminates 'Run' as a choice.\n")
    bob_available.pop(2) # Remove Run
    
    # --- Round 2 of Elimination ---
    print("Step 2: Analyzing the reduced 2x2 game (Rest vs. Bike).\n")
    
    # In the reduced game, Alice's payoffs are:
    # Rest: [0, 2]
    # Bike: [0, -2]
    # Rest (0, 2) weakly dominates Bike (0, -2) because 0==0 and 2>-2.
    # Alice will eliminate 'Bike'.
    print("For Alice: Comparing 'Rest' ([0, 2]) vs 'Bike' ([0, -2]) in the reduced game.")
    print("'Rest' payoffs are equal or better than 'Bike' payoffs.")
    print("=> Alice eliminates 'Bike' as a choice.\n")
    alice_available.pop(1) # Remove Bike

    # In the reduced game, Bob's payoffs are:
    # Rest: [0, 2]
    # Bike: [0, -2]
    # Rest (0, 2) weakly dominates Bike (0, -2).
    # Bob will eliminate 'Bike'.
    print("For Bob: Comparing 'Rest' ([0, 2]) vs 'Bike' ([0, -2]) in the reduced game.")
    print("'Rest' payoffs are equal or better than 'Bike' payoffs.")
    print("=> Bob eliminates 'Bike' as a choice.\n")
    bob_available.pop(1) # Remove Bike

    # Determine the final outcome
    final_alice_idx = list(alice_available.keys())[0]
    final_bob_idx = list(bob_available.keys())[0]

    final_alice_action = actions[final_alice_idx]
    final_bob_action = actions[final_bob_idx]
    final_payoff_tuple = payoffs[final_alice_idx][final_bob_idx]
    alice_utility = final_payoff_tuple[0]

    print("-" * 48)
    print("Conclusion:")
    print(f"After iterated elimination, the only rational outcome is for both to choose '{final_alice_action}'.")
    print(f"Predicted outcome: (Alice: {final_alice_action}, Bob: {final_bob_action})")
    
    print("\nCalculating Alice's Expected Utility:")
    print(f"Alice's utility is her payoff in the final (Rest, Rest) outcome.")
    print(f"Final Equation: Utility = Payoff_Alice(Rest, Rest) = {alice_utility}")
    print("\nAlice's final expected utility is:")
    print(alice_utility)
    
    # Final answer block as requested
    sys.stdout.write(f"\n<<<{alice_utility}>>>\n")

solve_utility()