import sys
# Redirect stdout to a string buffer to capture all prints
from io import StringIO
original_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """

    # Step 1 & 2: Define the payoffs for the symmetric outcomes, which are the only
    # ones considered by superrational agents.
    # If both rest, payoff is 0.
    # If both bike, payoff is -2.
    # If both run, payoff is -3.
    symmetric_payoffs = {
        "rest": 0,
        "bike": -2,
        "run": -3,
    }

    print("The problem assumes Alice and Bob are 'superrational'. This means they will each deduce the single best logical action and, being equally rational, will both choose the same action.")
    print("Therefore, we only need to consider the outcomes where they both make the same choice:\n")

    # Step 3: Evaluate and display the symmetric outcomes.
    rest_payoff = symmetric_payoffs["rest"]
    bike_payoff = symmetric_payoffs["bike"]
    run_payoff = symmetric_payoffs["run"]

    print(f"1. If both choose to rest, their payoff is: {rest_payoff}")
    print(f"2. If both choose to bike, their payoff is: {bike_payoff}")
    print(f"3. If both choose to run, their payoff is: {run_payoff}")

    # Step 4: A rational agent chooses the action that maximizes their payoff.
    # The superrational choice is the maximum of the symmetric payoffs.
    alice_utility = max(rest_payoff, bike_payoff, run_payoff)

    print("\nSince they are rational, they will choose the action that leads to the best possible outcome for themselves from this set.")
    
    # Step 5: Display the final calculation as an equation.
    # The problem asks to output the equation with each number.
    print("\nThe final calculation to determine this choice is finding the maximum of these payoffs:")
    final_equation = f"max({rest_payoff}, {bike_payoff}, {run_payoff}) = {alice_utility}"
    print(final_equation)

    print(f"\nThis means the logical choice for both Alice and Bob is to rest.")
    print(f"Therefore, Alice's expected utility is {alice_utility}.")

solve_utility()

# Get the captured output and print it to the actual console
output = captured_output.getvalue()
sys.stdout = original_stdout
print(output)
print("<<<0>>>")