import math

def solve_utility():
    """
    Calculates Alice's expected utility in a symmetric game with a superrational opponent.
    """
    # Define Alice's payoffs.
    # payoffs[alice_action][bob_action] = alice's_payoff
    payoffs = {
        'rest': {'rest': 0, 'bike': 2, 'run': 4},
        'bike': {'rest': 0, 'bike': -2, 'run': 2},
        'run':  {'rest': 0, 'bike': 0, 'run': -3}
    }

    actions = list(payoffs.keys())
    symmetric_payoffs = {}

    print("As a superrational agent in a symmetric game, Alice knows Bob will reason identically and make the same choice.")
    print("Therefore, she evaluates her payoff for each possible identical outcome:")

    for action in actions:
        # The payoff is determined by the diagonal of the game matrix.
        payoff = payoffs[action][action]
        symmetric_payoffs[action] = payoff
        print(f"- If both players choose to {action}, Alice's payoff is {payoff}.")

    # A rational agent maximizes their utility.
    # Find the maximum payoff among the symmetric outcomes.
    best_action = max(symmetric_payoffs, key=symmetric_payoffs.get)
    max_utility = symmetric_payoffs[best_action]

    payoff_values = list(symmetric_payoffs.values())
    
    # Construct the "equation" part of the output as requested.
    print("\nTo maximize her utility, Alice compares the payoffs of these outcomes.")
    print(f"The calculation is: max({payoff_values[0]}, {payoff_values[1]}, {payoff_values[2]}) = {max_utility}")

    print(f"\nThe best choice is to '{best_action}', resulting in a guaranteed utility of {max_utility}.")
    print("\nSince this outcome is certain under superrationality, Alice's expected utility is her utility in this single outcome.")
    print(f"Alice's expected utility is {max_utility}.")


solve_utility()
<<<0>>>