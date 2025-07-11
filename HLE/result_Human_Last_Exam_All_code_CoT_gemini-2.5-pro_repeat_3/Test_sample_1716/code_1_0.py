def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """

    # Define the payoff matrix from Alice's perspective.
    # The keys are Alice's actions, and the values are dictionaries
    # mapping Bob's action to Alice's payoff.
    payoffs = {
        "Rest": {"Rest": 0, "Bike": 2, "Run": 4},
        "Bike": {"Rest": 0, "Bike": -2, "Run": 2},
        "Run":  {"Rest": 0, "Bike": 0, "Run": -3}
    }
    
    actions = ["Rest", "Bike", "Run"]

    # In a superrational scenario, agents know they will independently
    # arrive at the same optimal strategy. Therefore, we only need to
    # consider the payoffs where both agents choose the same action.
    
    # Extract the payoffs from the diagonal of the matrix.
    diagonal_payoffs = {action: payoffs[action][action] for action in actions}
    
    print("Under superrationality, Alice and Bob know they will make the same choice.")
    print("We must find the best outcome for Alice among the scenarios where both act identically:")
    
    for action, payoff in diagonal_payoffs.items():
        print(f"- If both {action}, Alice's payoff is: {payoff}")
        
    # The optimal choice for a superrational agent is the one that maximizes
    # their payoff on the diagonal.
    max_payoff = max(diagonal_payoffs.values())
    
    # The expected utility is this maximum payoff, as the outcome is certain.
    
    # Create the string for the final equation
    payoff_values_str = ", ".join(map(str, diagonal_payoffs.values()))
    
    print("\nAlice's expected utility is the maximum of these possible payoffs.")
    print(f"Alice's Expected Utility = max({payoff_values_str})")
    print(f"Alice's Expected Utility = {max_payoff}")

solve_utility()