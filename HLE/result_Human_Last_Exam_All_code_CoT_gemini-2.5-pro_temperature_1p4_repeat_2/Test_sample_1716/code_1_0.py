def calculate_alices_utility():
    """
    Calculates Alice's expected utility in a superrational game.

    In a symmetric game with superrational agents, each agent knows the other
    will make the exact same choice. Therefore, we only need to evaluate the
    payoffs along the diagonal of the payoff matrix.
    """

    # Define Alice's payoffs for the three symmetric outcomes.
    # 1. Both Rest: Alice sees no one, payoff is 0.
    payoff_rest_rest = 0

    # 2. Both Bike: Alice gets angry at Bob for copying, payoff is -2.
    payoff_bike_bike = -2

    # 3. Both Run: Alice gets angry and suffers from the run, payoff is -3.
    payoff_run_run = -3

    # A superrational agent chooses the action that maximizes their payoff,
    # assuming the other agent will make the same choice.
    symmetric_payoffs = {
        "Rest": payoff_rest_rest,
        "Bike": payoff_bike_bike,
        "Run": payoff_run_run
    }

    # Find the action and utility from the best symmetric outcome.
    best_action = max(symmetric_payoffs, key=symmetric_payoffs.get)
    alices_expected_utility = symmetric_payoffs[best_action]

    print("This is a symmetric game between two superrational agents.")
    print("This means both Alice and Bob will independently choose the same action.")
    print("We must find which symmetric choice gives Alice the highest payoff:")
    print(f"\nIf both Rest, Alice's payoff is: {payoff_rest_rest}")
    print(f"If both Bike, Alice's payoff is: {payoff_bike_bike}")
    print(f"If both Run, Alice's payoff is: {payoff_run_run}")

    print("\nTo find her expected utility, Alice calculates the maximum of these potential outcomes.")
    print("The final equation for Alice's decision is:")
    # The final equation is printed here, showing each number.
    print(f"max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run}) = {alices_expected_utility}")

    print(f"\nThe best choice is to {best_action}.")
    print(f"Therefore, Alice's expected utility is {alices_expected_utility}.")


if __name__ == "__main__":
    calculate_alices_utility()