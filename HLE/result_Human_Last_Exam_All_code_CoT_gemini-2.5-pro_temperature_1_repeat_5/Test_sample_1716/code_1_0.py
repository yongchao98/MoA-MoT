def solve_superrational_game():
    """
    Calculates Alice's expected utility in a game with a superrational opponent.
    """
    # In a game with two identical, superrational players, they will both
    # deduce the same optimal strategy. This means we only need to consider
    # the outcomes where they both choose the same action.

    # 1. Define Alice's payoffs for the symmetric outcomes.
    payoff_rest = 0  # Payoff if both rest
    payoff_bike = -2 # Payoff if both bike
    payoff_run = -3  # Payoff if both run

    # 2. As a rational agent, Alice will choose the action that leads to the
    # maximum possible payoff among these symmetric outcomes.
    symmetric_payoffs = {
        "Rest": payoff_rest,
        "Bike": payoff_bike,
        "Run": payoff_run
    }

    best_choice = max(symmetric_payoffs, key=symmetric_payoffs.get)
    expected_utility = symmetric_payoffs[best_choice]

    # 3. Print the reasoning and the result.
    print("Alice and Bob are superrational, so they know they will make the same choice.")
    print("We only need to compare Alice's payoffs for the symmetric outcomes:")
    print(f"If both Rest, Alice's payoff is: {payoff_rest}")
    print(f"If both Bike, Alice's payoff is: {payoff_bike}")
    print(f"If both Run, Alice's payoff is: {payoff_run}")
    print("\nTo maximize her utility, Alice will choose the action with the highest payoff.")
    
    # The final equation showing the calculation.
    payoff_values = list(symmetric_payoffs.values())
    print(f"The calculation is: max({payoff_values[0]}, {payoff_values[1]}, {payoff_values[2]}) = {expected_utility}")
    
    print(f"\nTherefore, Alice's expected utility is {expected_utility}.")

solve_superrational_game()