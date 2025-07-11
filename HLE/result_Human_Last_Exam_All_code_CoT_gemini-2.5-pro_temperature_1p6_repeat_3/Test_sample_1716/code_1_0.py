import sys

def solve_superrational_game():
    """
    Calculates Alice's expected utility in a symmetric game between two
    superrational agents.
    """
    # Payoff matrix for Alice. The key is a tuple (Alice's Action, Bob's Action).
    # The value is Alice's resulting utility.
    payoffs = {
        ("Rest", "Rest"): 0,
        ("Bike", "Bike"): -2,
        ("Run", "Run"): -3,
        # Off-diagonal payoffs are not needed for superrational agents in a
        # symmetric game, but are included for completeness of the game definition.
        ("Rest", "Bike"): 2,
        ("Rest", "Run"): 4,
        ("Bike", "Rest"): 0,
        ("Bike", "Run"): 2,
        ("Run", "Rest"): 0,
        ("Run", "Bike"): 0
    }

    actions = ["Rest", "Bike", "Run"]

    print("Step 1: Identify the game type.")
    print("The game is symmetric, and the agents (Alice and Bob) are superrational.")
    print("This means they will both independently choose the same action.\n")

    print("Step 2: Identify the possible outcomes.")
    print("We only need to evaluate the outcomes where both agents pick the same action:")
    
    diagonal_payoffs = {action: payoffs[(action, action)] for action in actions}
    for action, payoff in diagonal_payoffs.items():
        print(f"- If both {action}, Alice's payoff is {payoff}")
    print("")

    print("Step 3: Determine the optimal choice.")
    print("Alice will choose the action that maximizes her payoff from this set of possibilities.")
    
    # Find the action with the maximum payoff
    best_action = max(diagonal_payoffs, key=diagonal_payoffs.get)
    max_utility = diagonal_payoffs[best_action]

    # Create the equation string as requested
    payoff_values = list(diagonal_payoffs.values())
    equation = f"max({payoff_values[0]}, {payoff_values[1]}, {payoff_values[2]}) = {max_utility}"
    
    print(f"The calculation is: {equation}\n")

    print("Step 4: State the conclusion.")
    print(f"The best choice for a superrational agent is '{best_action}', yielding a utility of {max_utility}.")
    print("Since this choice is deterministic, Alice's expected utility is this value.")

    # The final answer in the requested format is printed last
    # To avoid having the <<< >>>> in the middle of the text, we print it at the very end.
    # We use stderr for the explanation and stdout for the final answer block for clarity.
    # However, the user wants a single code block, so we'll just print it last.
    
    final_answer = f"<<<{max_utility}>>>"
    # This ensures the final answer is on the last line of the script output
    if 'site' in sys.modules:
        # In a typical execution environment
         print(f"\nFinal Answer: {final_answer}")
    else:
        # For platforms that might parse the output directly
        print(final_answer)


solve_superrational_game()