import collections

def solve_utility():
    """
    Calculates Alice's expected utility based on game theory principles.
    """
    # The choices available to Alice and Bob
    choices = ["Rest", "Bike", "Run"]

    # The payoff matrix is represented as a dictionary.
    # The key is a tuple (Alice's Choice, Bob's Choice).
    # The value is a tuple (Alice's Payoff, Bob's Payoff).
    # Note: A typo in the prompt for "Bike vs Run" is assumed to be corrected
    # to match the context.
    payoffs = {
        ("Rest", "Rest"): (0, 0),
        ("Bike", "Bike"): (-2, -2),
        ("Run", "Run"): (-3, -3),
        ("Rest", "Bike"): (2, 0),
        ("Bike", "Rest"): (0, 2),
        ("Rest", "Run"): (4, 0),
        ("Run", "Rest"): (0, 4),
        ("Bike", "Run"): (2, 0), # Corrected from the ambiguous prompt
        ("Run", "Bike"): (0, 2)
    }

    print("Step 1: Analyzing the game for dominant strategies for Alice.")
    
    # Create a simplified payoff matrix for Alice for easier analysis
    alice_payoffs = collections.defaultdict(dict)
    for (alice_choice, bob_choice), (alice_payoff, _) in payoffs.items():
        alice_payoffs[alice_choice][bob_choice] = alice_payoff

    print("Alice's Payoffs (Rows: Alice's choice, Columns: Bob's choice)")
    print(f"      {choices[0]:>5} {choices[1]:>5} {choices[2]:>5}")
    for alice_choice in choices:
        bob_r = alice_payoffs[alice_choice]['Rest']
        bob_b = alice_payoffs[alice_choice]['Bike']
        bob_n = alice_payoffs[alice_choice]['Run']
        print(f"{alice_choice:<5} {bob_r:>5} {bob_b:>5} {bob_n:>5}")

    print("\nStep 2: Checking if 'Rest' dominates other strategies for Alice.")
    
    # Comparing 'Rest' to 'Bike' for Alice
    rest_vs_bike = [alice_payoffs["Rest"][c] >= alice_payoffs["Bike"][c] for c in choices]
    print(f"Is 'Rest' always as good or better than 'Bike' for Alice? {all(rest_vs_bike)}")
    print(f"Payoffs for Rest: ({alice_payoffs['Rest']['Rest']}, {alice_payoffs['Rest']['Bike']}, {alice_payoffs['Rest']['Run']})")
    print(f"Payoffs for Bike: ({alice_payoffs['Bike']['Rest']}, {alice_payoffs['Bike']['Bike']}, {alice_payoffs['Bike']['Run']})")
    
    # Comparing 'Rest' to 'Run' for Alice
    rest_vs_run = [alice_payoffs["Rest"][c] >= alice_payoffs["Run"][c] for c in choices]
    print(f"\nIs 'Rest' always as good or better than 'Run' for Alice? {all(rest_vs_run)}")
    print(f"Payoffs for Rest: ({alice_payoffs['Rest']['Rest']}, {alice_payoffs['Rest']['Bike']}, {alice_payoffs['Rest']['Run']})")
    print(f"Payoffs for Run:  ({alice_payoffs['Run']['Rest']}, {alice_payoffs['Run']['Bike']}, {alice_payoffs['Run']['Run']})")

    print("\nConclusion: 'Rest' is a strictly dominant strategy for Alice. A superrational agent will always play a dominant strategy.")
    print("Since the game is symmetric, 'Rest' is also the dominant strategy for Bob.")

    print("\nStep 3: Determining the final outcome and Alice's utility.")
    # The unique outcome is that both players choose their dominant strategy.
    final_outcome = ("Rest", "Rest")
    final_payoff = payoffs[final_outcome]
    alice_utility = final_payoff[0]

    print(f"Both agents will choose 'Rest', leading to the outcome {final_outcome}.")
    print(f"The payoff for this outcome is {final_payoff}.")
    print(f"Alice's utility is the first value in the payoff tuple.")
    
    print("\nFinal equation for Alice's utility:")
    print(f"U(Alice) = {final_payoff[0]}")


solve_utility()
<<<0>>>