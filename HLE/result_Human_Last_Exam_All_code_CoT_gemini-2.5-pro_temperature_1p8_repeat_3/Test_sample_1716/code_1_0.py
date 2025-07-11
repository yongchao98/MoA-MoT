import collections

def solve_utility():
    """
    Calculates Alice's expected utility in a symmetric game with a superrational opponent.
    """
    # Step 1: Define the payoff matrix and actions
    # Actions are indexed: 0 for Rest, 1 for Bike, 2 for Run.
    # The matrix stores payoffs as (Alice's Utility, Bob's Utility).
    actions = ["Rest", "Bike", "Run"]
    payoffs = [
        # Bob's Actions: Rest,     Bike,      Run
        [(0, 0),   (2, 0),    (4, 0)],    # Alice's Action: Rest
        [(0, 2),   (-2, -2),  (2, 0)],    # Alice's Action: Bike
        [(0, 4),   (0, 2),    (-3, -3)],  # Alice's Action: Run
    ]

    print("Step 1: The payoff matrix for Alice and Bob is as follows:")
    print("        Bob: Rest   Bike      Run")
    for i, row in enumerate(payoffs):
        print(f"Alice: {actions[i]:<4} {str(row[0]):<9} {str(row[1]):<9} {str(row[2]):<9}")

    # Step 2: Apply the superrationality principle
    # Because the game is symmetric and both players are superrational, they must
    # independently choose the same action. We only need to examine the diagonal payoffs.
    print("\nStep 2: Apply the superrationality principle.")
    print("Alice knows Bob will reason identically and make the same choice.")
    print("Therefore, the only possible outcomes are when they both perform the same action:")
    
    alice_diagonal_payoffs = [payoffs[i][i][0] for i in range(len(actions))]
    
    for i, action in enumerate(actions):
        print(f"- If both players choose '{action}', Alice's utility is {alice_diagonal_payoffs[i]}.")

    # Step 3: Identify the optimal choice
    # Alice will choose the action that maximizes her utility from this set.
    print("\nStep 3: Alice chooses the action that maximizes her utility from these options.")
    print("We need to find the maximum value among Alice's possible outcomes.")
    
    # Step 4: Calculate the final utility by showing the equation
    # The problem asks to print the numbers in the final equation.
    print("\nFinal Calculation:")
    # Creates a string representation like "max(0, -2, -3)"
    equation_str = f"max({', '.join(map(str, alice_diagonal_payoffs))})"
    print(f"Alice's Expected Utility = {equation_str}")
    
    alice_expected_utility = max(alice_diagonal_payoffs)
    
    print(f"The result of this calculation is {alice_expected_utility}.")
    print(f"\nThus, Alice's expected utility is {alice_expected_utility}.")

solve_utility()