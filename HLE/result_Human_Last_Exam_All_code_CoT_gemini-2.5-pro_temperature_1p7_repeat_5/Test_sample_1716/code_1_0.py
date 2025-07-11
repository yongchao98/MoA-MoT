def solve_alice_utility():
    """
    Calculates Alice's expected utility based on the concept of superrationality.
    """
    # Step 1: Define the payoffs for the symmetric outcomes.
    # The problem states Alice and Bob are superrational agents. This means they know
    # that being identical reasoners, they will inevitably make the same choice.
    # Therefore, we only need to consider the scenarios where they both perform the same action.
    symmetric_payoffs = {
        'Rest': 0,  # Payoff if both rest
        'Bike': -2,  # Payoff if both bike
        'Run':  -3   # Payoff if both run
    }

    print("Alice and Bob are superrational, so they know they will make the same choice.")
    print("We analyze the payoffs for these symmetric outcomes:\n")

    for action, payoff in symmetric_payoffs.items():
        print(f"- If both choose to {action}, their payoff is {payoff}.")

    # Step 2: Find the action that leads to the best outcome among the symmetric choices.
    # The superrational choice is the one that maximizes their utility from this set.
    best_choice = max(symmetric_payoffs, key=symmetric_payoffs.get)
    alice_expected_utility = symmetric_payoffs[best_choice]

    print(f"\nThe optimal choice for both is '{best_choice}', as it yields the highest payoff.")

    # Step 3: Present the final calculation and result.
    # The final equation is finding the maximum of the possible symmetric payoffs.
    rest_payoff = symmetric_payoffs['Rest']
    bike_payoff = symmetric_payoffs['Bike']
    run_payoff = symmetric_payoffs['Run']

    print("\nThe final calculation is:")
    print(f"Alice's Expected Utility = max({rest_payoff}, {bike_payoff}, {run_payoff})")
    print(f"Alice's expected utility is {alice_expected_utility}.")


solve_alice_utility()