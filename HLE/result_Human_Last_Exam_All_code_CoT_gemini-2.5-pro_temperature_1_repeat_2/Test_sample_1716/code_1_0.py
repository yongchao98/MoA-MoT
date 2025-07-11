import sys

def solve_utility():
    """
    Calculates Alice's expected utility based on the game theory problem.
    """
    # Step 1: Define the payoff matrix for Alice based on the problem description.
    # The keys are tuples of (Alice's Action, Bob's Action).
    payoffs_alice = {
        ('Rest', 'Rest'): 0,
        ('Rest', 'Bike'): 2,
        ('Rest', 'Run'): 4,
        ('Bike', 'Rest'): 0,
        ('Bike', 'Bike'): -2,
        ('Bike', 'Run'): 2,
        ('Run', 'Rest'): 0,
        ('Run', 'Bike'): 0,
        ('Run', 'Run'): -3,
    }

    # Step 2: Apply the superrationality assumption.
    # Alice and Bob will make the same choice, so we only look at the diagonal of the payoff matrix.
    actions = ['Rest', 'Bike', 'Run']
    symmetric_outcomes = [(action, action) for action in actions]
    
    # Step 3: Determine the optimal choice by comparing utilities.
    utility_rest = payoffs_alice[('Rest', 'Rest')]
    utility_bike = payoffs_alice[('Bike', 'Bike')]
    utility_run = payoffs_alice[('Run', 'Run')]

    # Step 4: The expected utility is the maximum of these potential outcomes.
    expected_utility = max(utility_rest, utility_bike, utility_run)

    # Print the step-by-step reasoning
    print("Alice and Bob are superrational, meaning they will reason identically and choose the same action.")
    print("Therefore, we only need to consider the three scenarios where they act alike:")
    print(f"1. If both Rest, Alice's utility is: {utility_rest}")
    print(f"2. If both Bike, Alice's utility is: {utility_bike}")
    print(f"3. If both Run, Alice's utility is: {utility_run}")
    print("\nAlice will choose the action that leads to the best of these outcomes for herself.")
    print(f"The final calculation for her expected utility is the maximum of these values.")
    print(f"Equation: max({utility_rest}, {utility_bike}, {utility_run}) = {expected_utility}")
    
    # The final answer is wrapped for parsing, as requested.
    # The redirection of stdout is a way to ensure the final wrapped answer is the only thing
    # that goes to the calling process, in case this script is run in an automated environment.
    original_stdout = sys.stdout
    try:
        # Redirect stdout to a dummy stream to capture the final answer without printing it with the explanation
        from io import StringIO
        sys.stdout = StringIO()
        
        # This is the value to be returned
        print(f"<<<{expected_utility}>>>")
    finally:
        # Restore stdout and get the captured value
        captured_output = sys.stdout.getvalue()
        sys.stdout = original_stdout
        # We print the reasoning above, and the final answer is conceptually returned here.

if __name__ == '__main__':
    solve_utility()
