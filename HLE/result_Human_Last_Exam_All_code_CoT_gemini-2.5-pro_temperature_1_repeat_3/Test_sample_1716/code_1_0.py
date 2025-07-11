import sys

def solve_superrational_game():
    """
    Solves for Alice's expected utility in a superrational game.
    """
    # Step 1: Define the payoff matrix for Alice.
    # We only need Alice's payoffs since Bob's reasoning will mirror hers.
    # The matrix is represented as a dictionary: Alice's Action -> {Bob's Action: Alice's Payoff}
    payoffs_alice = {
        'Rest': {'Rest': 0, 'Bike': 2, 'Run': 4},
        'Bike': {'Rest': 0, 'Bike': -2, 'Run': 2},
        'Run':  {'Rest': 0, 'Bike': 0, 'Run': -3}
    }
    
    actions = ['Rest', 'Bike', 'Run']

    # Step 2: Explain and apply the concept of superrationality.
    print("Superrationality implies that Alice knows Bob will make the exact same decision as her.")
    print("Therefore, we only need to consider the outcomes where both choose the same action:")
    
    # Step 3: Determine the optimal choice from the symmetric outcomes.
    symmetric_outcomes = {}
    for action in actions:
        payoff = payoffs_alice[action][action]
        symmetric_outcomes[action] = payoff
        print(f"- If both Rest: Alice's payoff is {payoffs_alice['Rest']['Rest']}")
        print(f"- If both Bike: Alice's payoff is {payoffs_alice['Bike']['Bike']}")
        print(f"- If both Run:  Alice's payoff is {payoffs_alice['Run']['Run']}")
        break # Only need to print this explanation once.

    # Find the maximum possible payoff for Alice under this constraint.
    best_payoff = -sys.maxsize
    best_action = None
    
    diagonal_payoffs = []
    for action in actions:
        payoff = symmetric_outcomes[action]
        diagonal_payoffs.append(payoff)
        if payoff > best_payoff:
            best_payoff = payoff
            best_action = action
            
    # Step 4: Print the final calculation and the resulting expected utility.
    print("\nAlice must choose the action that maximizes her payoff, assuming Bob does the same.")
    print("The final equation is finding the maximum of these potential outcomes:")
    
    # Format the numbers into the equation string as requested.
    payoff_str_list = [str(p) for p in diagonal_payoffs]
    equation = f"max({', '.join(payoff_str_list)}) = {best_payoff}"
    
    print(equation)
    
    print(f"\nThe best choice for both is '{best_action}', leading to a payoff of {best_payoff}.")
    print(f"Thus, Alice's expected utility is {best_payoff}.")

solve_superrational_game()
print("<<<0>>>")