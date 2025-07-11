def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """
    # Step 1: Define Alice's payoff matrix.
    # The matrix is structured as:
    #           Bob: Rest | Bob: Bike | Bob: Run
    # Alice: Rest      ...       ...        ...
    # Alice: Bike      ...       ...        ...
    # Alice: Run       ...       ...        ...
    payoffs_alice = {
        "Rest": {"Rest": 0, "Bike": 2, "Run": 4},
        "Bike": {"Rest": 0, "Bike": -2, "Run": 2},
        "Run":  {"Rest": 0, "Bike": 0, "Run": -3}
    }

    # Step 2: Apply the superrationality constraint.
    # Alice knows that Bob, being her superrational counterpart, will make the same choice as her.
    # Therefore, we only need to consider the diagonal of the payoff matrix.
    print("Based on the principle of superrationality, both Alice and Bob will choose the same action.")
    print("We must therefore evaluate the payoffs for the following scenarios:")

    payoff_rest_rest = payoffs_alice["Rest"]["Rest"]
    payoff_bike_bike = payoffs_alice["Bike"]["Bike"]
    payoff_run_run = payoffs_alice["Run"]["Run"]

    print(f"- If both Rest, Alice's payoff is: {payoff_rest_rest}")
    print(f"- If both Bike, Alice's payoff is: {payoff_bike_bike}")
    print(f"- If both Run, Alice's payoff is: {payoff_run_run}")
    print("\n")

    # Step 3 & 4: Identify the optimal choice and determine the expected utility.
    # Alice will choose the action that maximizes her payoff from this set of possibilities.
    # Her expected utility is the payoff from that single, deterministic outcome.
    expected_utility = max(payoff_rest_rest, payoff_bike_bike, payoff_run_run)

    print("Alice will choose the action that results in the best outcome for herself from these options.")
    print("Alice's expected utility is the maximum of these potential payoffs.")
    print(f"The final calculation is: max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run}) = {expected_utility}")
    print("\n")
    print(f"Therefore, Alice's expected utility is {expected_utility}.")


solve_utility()