def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """
    # The actions are ordered as: Rest, Bike, Run
    actions = ["Rest", "Bike", "Run"]

    # The payoff matrix is defined as a list of lists.
    # Each cell is a tuple: (Alice's Payoff, Bob's Payoff)
    # Alice's choices are rows, Bob's choices are columns.
    payoff_matrix = [
        # Bob: Rest    Bike        Run
        [(0, 0),    (2, 0),    (4, 0)],   # Alice: Rest
        [(0, 2),    (-2, -2),  (2, 0)],   # Alice: Bike
        [(0, 4),    (0, 2),    (-3, -3)]   # Alice: Run
    ]

    print("Step 1: Analyze symmetric outcomes due to superrationality.")
    print("Since Alice and Bob are superrational, they know they will reason identically and choose the same action.")
    print("Therefore, we only need to consider the diagonal of the payoff matrix:")

    symmetric_payoffs = []
    for i in range(len(actions)):
        # Get Alice's payoff when both choose the same action
        payoff = payoff_matrix[i][i][0]
        symmetric_payoffs.append(payoff)
        print(f"If both players choose '{actions[i]}', Alice's payoff is: {payoff}")

    print("\nStep 2: Determine the optimal choice.")
    print("The agents will choose the action that leads to the best of these symmetric outcomes.")

    # To satisfy the request "output each number in the final equation"
    # we format the calculation into a string.
    equation_str = "Alice's expected utility = max(" + ", ".join(map(str, symmetric_payoffs)) + ")"
    print(equation_str)

    best_utility = max(symmetric_payoffs)

    print(f"\nStep 3: State the final utility.")
    print(f"The maximum payoff Alice can expect in a symmetric outcome is {best_utility}.")

solve_utility()
<<<0>>>