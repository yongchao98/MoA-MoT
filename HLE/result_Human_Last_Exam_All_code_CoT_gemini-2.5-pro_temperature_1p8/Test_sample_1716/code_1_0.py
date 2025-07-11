import numpy as np

def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """

    # Define the choices available to Alice and Bob
    choices = ["Rest", "Bike", "Run"]

    # Define Alice's payoff matrix based on the problem description.
    # Rows represent Alice's choice, Columns represent Bob's choice.
    # payoffs_alice[i, j] is Alice's payoff if she chooses i and Bob chooses j.
    payoffs_alice = np.array([
      # Bob: Rest, Bike, Run
        [0,      2,      4],    # Alice: Rest
        [0,     -2,      2],    # Alice: Bike
        [0,      0,     -3]     # Alice: Run
    ])

    print("--- Step 1: Understanding the Payoff Matrix ---")
    print("Alice's payoffs are structured as follows (Rows: Alice's choice, Columns: Bob's choice):")
    print("Choices: {}".format(choices))
    print(payoffs_alice)
    print("\n")

    print("--- Step 2: Applying Superrationality ---")
    print("Alice and Bob are 'superrational' and 'extremely self-aware'. This means they know they will both independently reason to the same conclusion.")
    print("Therefore, Alice only needs to consider the scenarios where she and Bob make the same choice.")
    print("These scenarios correspond to the diagonal of the payoff matrix:\n")

    # Extract the payoffs from the diagonal of the matrix
    symmetric_payoffs = np.diag(payoffs_alice)
    
    payoff_rest_rest = payoffs_alice[0, 0]
    payoff_bike_bike = payoffs_alice[1, 1]
    payoff_run_run = payoffs_alice[2, 2]

    print(f"1. If both Rest, Alice's payoff is {payoff_rest_rest}")
    print(f"2. If both Bike, Alice's payoff is {payoff_bike_bike}")
    print(f"3. If both Run, Alice's payoff is {payoff_run_run}")
    print("\n")
    
    print("--- Step 3: Calculating Expected Utility ---")
    print("To maximize her utility, Alice will choose the action that results in the best payoff from these symmetric outcomes.")
    print("This is equivalent to finding the maximum value among the payoffs for these choices.")

    # Calculate the maximum of the symmetric payoffs
    expected_utility = np.max(symmetric_payoffs)

    # Print the final equation as requested
    print("\nFinal Calculation:")
    print(f"Alice's Expected Utility = max(Payoff(Rest,Rest), Payoff(Bike,Bike), Payoff(Run,Run))")
    print(f"Alice's Expected Utility = max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run}) = {expected_utility}")
    print("\n")
    print(f"Thus, Alice's expected utility is {expected_utility}.")


solve_utility()