def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """

    # According to the problem description, these are the payoffs for Alice
    # if both agents choose the same action.
    payoff_if_both_rest = 0
    payoff_if_both_bike = -2
    payoff_if_both_run = -3

    # In a game with superrational agents, each knows the other will use the
    # same logic. Therefore, they will arrive at the same action. We only
    # need to consider these symmetric outcomes.
    symmetric_payoffs = {
        "rest": payoff_if_both_rest,
        "bike": payoff_if_both_bike,
        "run": payoff_if_both_run
    }

    # Alice, being rational, will choose the action that maximizes her utility
    # from this set of possibilities.
    alices_choice = max(symmetric_payoffs, key=symmetric_payoffs.get)
    expected_utility = symmetric_payoffs[alices_choice]

    print("Alice and Bob are superrational, so they will both choose the same action.")
    print("We must find which symmetric outcome is best for Alice:")
    print(f"- If both Rest, Alice's utility is {payoff_if_both_rest}.")
    print(f"- If both Bike, Alice's utility is {payoff_if_both_bike}.")
    print(f"- If both Run, Alice's utility is {payoff_if_both_run}.")
    print("\nTo find her expected utility, Alice calculates the maximum of these potential outcomes.")
    # The final equation demonstrates the comparison of the symmetric payoffs.
    print(f"Final calculation: max({payoff_if_both_rest}, {payoff_if_both_bike}, {payoff_if_both_run}) = {expected_utility}")
    print(f"\nAlice's expected utility is {expected_utility}.")

solve_utility()
<<<0>>>