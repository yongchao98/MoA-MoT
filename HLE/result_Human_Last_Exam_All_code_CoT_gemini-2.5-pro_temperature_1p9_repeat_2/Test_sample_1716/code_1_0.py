import math

def solve_utility():
    """
    Solves for Alice's expected utility in a superrational game.
    """
    # 1. Define the actions and the payoff matrix based on the problem description.
    # The matrix stores tuples of (Alice's payoff, Bob's payoff).
    payoff_matrix = {
        #           Bob's choice
        # Alice's  |   rest   |   bike   |    run
        #----------------------------------------------
        'rest':    { 'rest':(0,0),  'bike':(2,0),  'run':(4,0) },
        'bike':    { 'rest':(0,2),  'bike':(-2,-2), 'run':(2,0) },
        'run':     { 'rest':(0,4),  'bike':(0,2),  'run':(-3,-3)}
    }

    actions = ['rest', 'bike', 'run']

    # 2. Apply the principle of superrationality.
    # Superrational players know they will make the same choice.
    # We only need to examine the outcomes where their actions are identical.
    symmetric_payoffs = {}
    for action in actions:
        # Get Alice's payoff for each symmetric outcome
        symmetric_payoffs[action] = payoff_matrix[action][action][0]

    # 3. Determine the optimal choice for a superrational player.
    # They will choose the action that leads to the highest payoff among the symmetric options.
    max_utility = -math.inf
    best_action = None

    for action, payoff in symmetric_payoffs.items():
        if payoff > max_utility:
            max_utility = payoff
            best_action = action
            
    # 4. Print the reasoning and the final answer.
    print("Step 1: Analyzing the payoffs for symmetric outcomes (where Alice and Bob choose the same action):")
    for action, payoff in symmetric_payoffs.items():
        print(f"  - If both {action}, Alice's payoff is: {payoff}")
    
    print("\nStep 2: A superrational player knows the other will make the same choice, so they will choose the action that maximizes their own payoff from the symmetric outcomes.")
    
    payoff_values = list(symmetric_payoffs.values())
    print(f"\nStep 3: Calculating the maximum utility from the set of possible payoffs {payoff_values}.")
    
    # The 'equation' is finding the maximum of the symmetric payoffs.
    equation_str = "max(" + ", ".join(map(str, payoff_values)) + ")"
    
    print(f"Final Equation: Alice's Expected Utility = {equation_str} = {max_utility}")


solve_utility()
<<<0>>>