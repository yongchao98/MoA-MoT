import sys

def solve_utility():
    """
    Calculates Alice's expected utility based on superrationality in a symmetric game.
    """
    # Step 1: Define the payoffs for Alice in the scenarios where both players choose the same action.
    # According to the problem description:
    # - If both rest, payoff is 0.
    # - If both bike, payoff is -2.
    # - If both run, payoff is -3.
    payoff_rest_rest = 0
    payoff_bike_bike = -2
    payoff_run_run = -3

    # Step 2: Explain the reasoning based on superrationality.
    # Alice knows Bob is also superrational and is in a symmetric situation.
    # This means Bob will use the same logic and arrive at the same action as Alice.
    # Therefore, they only need to consider the outcomes where they both choose the same action.
    print("Thinking Process:")
    print("1. Alice and Bob are superrational agents in a symmetric game.")
    print("2. They will both deduce the same optimal strategy, so they will choose the same action.")
    print("3. They must evaluate the outcomes where their actions are identical:")
    print(f"   - If both Rest, Alice's payoff is {payoff_rest_rest}.")
    print(f"   - If both Bike, Alice's payoff is {payoff_bike_bike}.")
    print(f"   - If both Run, Alice's payoff is {payoff_run_run}.")
    print("4. Alice will choose the action that leads to the best payoff for herself from this set.")
    print("\nCalculating the result...")

    # Step 3: Find the optimal choice by selecting the maximum of these payoffs.
    # This represents the outcome both agents will choose.
    alice_utility = max(payoff_rest_rest, payoff_bike_bike, payoff_run_run)

    # Step 4: Print the final calculation and the answer.
    # The final equation shows the comparison being made.
    print("The final utility is the maximum of the possible symmetric outcomes:")
    print(f"max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run}) = {alice_utility}")
    print(f"\nThus, both Alice and Bob will choose to Rest.")
    print(f"Alice's expected utility is {alice_utility}.")

    # Output the final answer in the required format
    sys.stdout.write(f"\n<<<{float(alice_utility)}>>>")

solve_utility()