import sys

def solve_superrational_game():
    """
    Calculates Alice's expected utility in a superrational game.
    """
    # Step 1 & 2: Define the payoffs for Alice when she and Bob choose the same action.
    # These are the only outcomes possible under the superrationality assumption.
    payoff_both_rest = 0
    payoff_both_bike = -2
    payoff_both_run = -3

    # Create a dictionary to hold the payoffs for clarity.
    payoffs = {
        "Rest": payoff_both_rest,
        "Bike": payoff_both_bike,
        "Run": payoff_both_run
    }

    print("Analyzing the game for a superrational agent...")
    print("Superrational agents know their counterpart will reason identically and make the same choice.")
    print("Therefore, we only need to consider the payoffs where both Alice and Bob choose the same action:\n")

    for action, payoff in payoffs.items():
        print(f"  - If both choose to {action}, Alice's payoff is: {payoff}")

    # Step 3: A rational agent chooses the action that maximizes their utility.
    # We find the maximum payoff among the possible outcomes.
    # The .get method on the dictionary allows us to find the key (action) for the max value.
    optimal_choice = max(payoffs, key=payoffs.get)
    expected_utility = payoffs[optimal_choice]

    # Step 4: Print the reasoning and the final calculation for expected utility.
    print("\nTo maximize her utility, Alice will choose the action with the highest payoff.")
    print("The calculation is to find the maximum of the possible payoffs.")

    # This fulfills the requirement to output each number in the final equation.
    payoff_values = list(payoffs.values())
    print(f"Equation: max({payoff_values[0]}, {payoff_values[1]}, {payoff_values[2]}) = {expected_utility}")

    print(f"\nThe optimal choice for both Alice and Bob is to {optimal_choice}.")
    print(f"Therefore, Alice's expected utility is {expected_utility}.")

if __name__ == "__main__":
    solve_superrational_game()