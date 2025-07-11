def solve_utility():
    """
    Calculates Alice's expected utility in a game with a superrational opponent.
    """
    # Step 1: Define the payoff matrix based on the problem description.
    # The matrix is represented as a dictionary where keys are (Alice's Action, Bob's Action)
    # and values are (Alice's Payoff, Bob's Payoff).
    payoffs = {
        ('Rest', 'Rest'): (0, 0),
        ('Rest', 'Bike'): (2, 0),
        ('Rest', 'Run'): (4, 0),
        ('Bike', 'Rest'): (0, 2),
        ('Bike', 'Bike'): (-2, -2),
        ('Bike', 'Run'): (2, 0),
        ('Run', 'Rest'): (0, 4),
        ('Run', 'Bike'): (0, 2),
        ('Run', 'Run'): (-3, -3)
    }
    actions = ['Rest', 'Bike', 'Run']

    # Step 2 & 3: Apply the principle of superrationality.
    # Superrational agents know they will reason identically and thus choose the same action.
    # We only need to consider the symmetric outcomes on the diagonal of the payoff matrix.
    print("Alice and Bob are superrational, so they know they will both choose the same action.")
    print("Therefore, they only need to consider the outcomes where their actions are identical:")

    symmetric_payoffs = {}
    for action in actions:
        payoff = payoffs[(action, action)][0] # Alice's payoff
        symmetric_payoffs[action] = payoff
        print(f"- If both choose to {action}, Alice's payoff is {payoff}.")

    # Step 4: Determine the optimal choice by finding the best symmetric outcome.
    # Alice will choose the action that leads to the highest payoff for herself from this set.
    best_action = max(symmetric_payoffs, key=symmetric_payoffs.get)
    max_payoff = symmetric_payoffs[best_action]

    # Step 5: The expected utility is the payoff from the single, determined outcome.
    print("\nTo maximize her utility, Alice must choose the action with the highest payoff among these options.")
    
    # The final equation is a comparison of the possible payoffs.
    payoff_values = list(symmetric_payoffs.values())
    print(f"The final decision is based on comparing the payoffs: {payoff_values[0]}, {payoff_values[1]}, and {payoff_values[2]}.")
    print(f"The maximum value is {max_payoff}, which corresponds to the action '{best_action}'.")
    
    print(f"\nSince both agents will choose to '{best_action}', this outcome is certain.")
    print(f"Alice's expected utility is her payoff in this outcome.")
    print(f"Alice's Expected Utility = {max_payoff}")

solve_utility()