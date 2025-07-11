import numpy as np

def solve_utility():
    """
    Solves the game theory problem for Alice and Bob.
    """
    # Step 1: Define the payoff matrix (Alice's payoff, Bob's payoff)
    # The actions are ordered as 'Rest', 'Bike', 'Run' for both rows (Alice) and columns (Bob)
    payoffs = {
        'Rest': {'Rest': (0, 0), 'Bike': (2, 0), 'Run': (4, 0)},
        'Bike': {'Rest': (0, 2), 'Bike': (-2, -2), 'Run': (2, 0)},
        'Run':  {'Rest': (0, 4), 'Bike': (0, 2), 'Run': (-3, -3)}
    }
    
    actions = ['Rest', 'Bike', 'Run']
    
    # Extract Alice's payoffs into a more accessible structure
    alice_payoffs = {action: [] for action in actions}
    for alice_action in actions:
        for bob_action in actions:
            alice_payoffs[alice_action].append(payoffs[alice_action][bob_action][0])

    print("Analyzing the game from Alice's perspective...")
    print("Alice's payoffs if Bob chooses Rest, Bike, or Run respectively:")
    for action in actions:
        print(f"  - If Alice chooses {action}: {alice_payoffs[action]}")
    print("-" * 30)

    # Step 2 & 3: Find weakly dominated strategies for Alice
    print("Checking for weakly dominated strategies...")
    
    dominated_strategies = set()
    
    # Compare every pair of strategies to find dominance
    for i, action1 in enumerate(actions):
        for j, action2 in enumerate(actions):
            if i == j:
                continue
            
            payoffs1 = np.array(alice_payoffs[action1])
            payoffs2 = np.array(alice_payoffs[action2])
            
            # Check if payoffs1 is always >= payoffs2
            if np.all(payoffs1 >= payoffs2) and np.any(payoffs1 > payoffs2):
                print(f"Alice's strategy '{action1}' weakly dominates '{action2}'.")
                dominated_strategies.add(action2)

    if not dominated_strategies:
        print("No weakly dominated strategies found.")
        return

    # Step 4: Determine the rational choice by eliminating dominated strategies
    rational_choices = [action for action in actions if action not in dominated_strategies]
    
    print("\nAs superrational agents, both Alice and Bob will eliminate their weakly dominated strategies.")
    print(f"The eliminated strategies are: {list(dominated_strategies)}")
    
    if len(rational_choices) == 1:
        final_choice = rational_choices[0]
        print(f"This leaves '{final_choice}' as the only rational choice for both.")
        
        # Step 5: Calculate Alice's final utility
        # The outcome is (final_choice, final_choice) because Bob reasons symmetrically
        final_outcome_payoff = payoffs[final_choice][final_choice]
        alices_utility = final_outcome_payoff[0]
        
        print("\nThe predicted outcome is (Alice: {}, Bob: {}).".format(final_choice, final_choice))
        print("The payoff for this outcome is {}.".format(final_outcome_payoff))
        
        # The final equation is simply the value from the payoff matrix for the determined outcome.
        print(f"\nAlice's expected utility = Payoff_Alice({final_choice}, {final_choice}) = {alices_utility}")
        print(f"<<<{float(alices_utility)}>>>")
    else:
        print("Analysis is inconclusive as multiple rational choices remain.")


solve_utility()