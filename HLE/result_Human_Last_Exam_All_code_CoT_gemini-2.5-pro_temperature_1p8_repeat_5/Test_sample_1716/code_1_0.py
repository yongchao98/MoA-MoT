import collections

def solve_superrational_game():
    """
    Solves the game theory problem for Alice's expected utility
    assuming superrational agents.
    """

    # Step 1: Model the game with a payoff matrix for Alice.
    # The first key is Alice's action, the second key is Bob's action.
    alice_payoffs = collections.OrderedDict([
        ('Rest', {'Rest': 0, 'Bike': 2, 'Run': 4}),
        ('Bike', {'Rest': 0, 'Bike': -2, 'Run': 2}),
        ('Run',  {'Rest': 0, 'Bike': 0, 'Run': -3})
    ])
    actions = list(alice_payoffs.keys())

    # Step 2: Explain the superrationality principle.
    print("Analyzing the game for superrational agents Alice and Bob.")
    print("A superrational agent knows that an equally superrational opponent will use the same logic and arrive at the same conclusion.")
    print("Therefore, both Alice and Bob will choose the same action.")
    print("\nWe only need to evaluate the outcomes on the 'diagonal' of the payoff matrix:")

    # Step 3: Identify the payoffs for symmetric choices.
    symmetric_payoffs = collections.OrderedDict()
    for action in actions:
        payoff = alice_payoffs[action][action]
        symmetric_payoffs[action] = payoff
        print(f"If both players {action}, Alice's payoff is {payoff}.")

    # Step 4: Find the optimal choice that maximizes utility.
    max_payoff = -float('inf')
    best_action = None

    for action, payoff in symmetric_payoffs.items():
        if payoff > max_payoff:
            max_payoff = payoff
            best_action = action

    # Step 5: Print the final calculation and conclusion.
    print("\nTo maximize her utility, Alice will choose the action that gives her the best payoff from these symmetric outcomes.")
    
    payoff_numbers = list(symmetric_payoffs.values())
    print(f"The equation is to find the maximum of these payoffs: max({payoff_numbers[0]}, {payoff_numbers[1]}, {payoff_numbers[2]})")
    
    print(f"The maximum payoff is {max_payoff}, which comes from the action '{best_action}'.")
    print("\nThus, Alice will choose to Rest, and she knows Bob will do the same.")
    print(f"\nAlice's expected utility is: {max_payoff}")

solve_superrational_game()
<<<0>>>