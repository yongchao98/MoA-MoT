import numpy as np

def solve_utility():
    """
    Calculates Alice's expected utility based on the principles of superrationality.
    """
    # Step 1: Define the payoff matrix for Alice based on the problem description.
    # Rows correspond to Alice's choice (0:Rest, 1:Bike, 2:Run)
    # Columns correspond to Bob's choice (0:Rest, 1:Bike, 2:Run)
    payoffs_alice = np.array([
        # Bob: Rest, Bike, Run
        [0,  2,  4],  # Alice: Rest
        [0, -2,  2],  # Alice: Bike
        [0,  0, -3]   # Alice: Run
    ])

    print("Step 1: The payoff matrix for Alice is:")
    print(payoffs_alice)
    print("\n-------------------------------------\n")

    # Step 2: Apply the principle of superrationality.
    # Because Alice and Bob are superrational and the game is symmetric,
    # they will both independently arrive at the same decision.
    # This means we only need to consider the outcomes where they both pick the same action.
    # These are the payoffs on the main diagonal of the matrix.
    symmetric_payoffs = np.diag(payoffs_alice)
    actions = ["Rest", "Bike", "Run"]

    print("Step 2: Under superrationality, only symmetric outcomes are possible.")
    print(f"If both Rest, Alice's utility is: {symmetric_payoffs[0]}")
    print(f"If both Bike, Alice's utility is: {symmetric_payoffs[1]}")
    print(f"If both Run, Alice's utility is: {symmetric_payoffs[2]}")
    print("\n-------------------------------------\n")


    # Step 3: Alice chooses the action that maximizes her utility among these options.
    # This is her expected utility.
    expected_utility = np.max(symmetric_payoffs)

    print("Step 3: Alice will choose the action that leads to the best outcome for her.")
    print("This is the maximum of the possible symmetric payoffs.")
    print("\nFinal Calculation:")
    # The final print statement shows the equation with all numbers.
    print(f"Alice's expected utility = max({symmetric_payoffs[0]}, {symmetric_payoffs[1]}, {symmetric_payoffs[2]}) = {expected_utility}")

solve_utility()