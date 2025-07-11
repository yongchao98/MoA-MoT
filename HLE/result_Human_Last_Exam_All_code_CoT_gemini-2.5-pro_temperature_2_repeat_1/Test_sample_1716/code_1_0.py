import sys

def solve_utility():
    """
    Calculates Alice's expected utility in a game with a superrational opponent.
    """
    # Step 1: Establish the payoff matrix.
    # The matrix is represented as a list of lists.
    # Each inner list is a row for Alice's choice (Rest, Bike, Run).
    # Each tuple is the payoff for (Alice, Bob).
    payoff_matrix = [
        # Bob's choice: Rest,      Bike,       Run
        [(0, 0),    (2, 0),     (4, 0)],    # Alice rests
        [(0, 2),    (-2, -2),   (2, 0)],    # Alice bikes
        [(0, 4),    (0, 2),     (-3, -3)]    # Alice runs
    ]

    # Step 2 & 3: Apply the "superrational" concept.
    # Superrational players in a symmetric game will independently arrive at the
    # same strategy. We only need to evaluate the outcomes on the main diagonal
    # where both players make the same choice.

    # Alice's payoff if both players Rest (row 0, col 0)
    payoff_both_rest = payoff_matrix[0][0][0]

    # Alice's payoff if both players Bike (row 1, col 1)
    payoff_both_bike = payoff_matrix[1][1][0]

    # Alice's payoff if both players Run (row 2, col 2)
    payoff_both_run = payoff_matrix[2][2][0]
    
    # Step 4: Determine the optimal choice by comparing these payoffs.
    # A rational agent will choose the action that maximizes their payoff.
    alice_utility = max(payoff_both_rest, payoff_both_bike, payoff_both_run)

    # Step 5 & 6: Print the reasoning and the result.
    print("Based on the superrationality assumption, Alice knows Bob will reason identically and make the same choice.")
    print("Therefore, we compare the payoffs for Alice from the symmetric outcomes:")
    print(f"1. If both Rest, Alice's utility is: {payoff_both_rest}")
    print(f"2. If both Bike, Alice's utility is: {payoff_both_bike}")
    print(f"3. If both Run, Alice's utility is: {payoff_both_run}")
    print("\nThe optimal choice is the one with the maximum utility.")
    print(f"The final calculation is finding the maximum of these values:")
    print(f"max({payoff_both_rest}, {payoff_both_bike}, {payoff_both_run}) = {alice_utility}")
    print(f"\nThus, both will choose to Rest, and Alice's expected utility is {alice_utility}.")

    # Appending the final answer as per instructions
    sys.stdout.write("\n<<<")
    sys.stdout.write(str(alice_utility))
    sys.stdout.write(">>>\n")

solve_utility()