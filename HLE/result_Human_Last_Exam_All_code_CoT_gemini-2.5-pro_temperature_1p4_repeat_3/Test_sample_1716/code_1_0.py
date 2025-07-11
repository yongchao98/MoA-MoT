import numpy as np

def solve_utility():
    """
    Calculates Alice's expected utility in a game with a superrational opponent.
    """
    actions = ["Rest", "Bike", "Run"]
    
    # Payoff matrix for Alice. U[i][j] is Alice's payoff if she chooses action i and Bob chooses action j.
    # Rows: Alice's actions (Rest, Bike, Run)
    # Columns: Bob's actions (Rest, Bike, Run)
    # Payoffs derived from the problem description:
    # (R,R)=(0,0), (B,B)=(-2,-2), (N,N)=(-3,-3)
    # (R,B)=(2,0), (R,N)=(4,0)
    # (B,R)=(0,2), (B,N)=(2,0)
    # (N,R)=(0,4), (N,B)=(0,2)
    # Note: The problem had a typo for (B,R) which was inferred to be for (B,N)
    utility_matrix_alice = np.array([
        [0, 2, 4],  # Alice Rests
        [0, -2, 2],  # Alice Bikes
        [0, 0, -3]   # Alice Runs
    ])

    print("Analyzing the game for superrational players Alice and Bob.")
    print("The key insight is that superrational players in a symmetric game will play the same strategy.")
    print("This means we must find the best symmetric Nash Equilibrium.\n")
    print("Alice's Payoff Matrix (Rows=Alice's choice, Cols=Bob's choice):")
    print(utility_matrix_alice)
    print("\nWe will test each of the three pure symmetric strategies (Rest,Rest), (Bike,Bike), (Run,Run) to see if they are Nash Equilibria.\n")

    symmetric_equilibria = []

    # Check each pure strategy on the diagonal of the payoff matrix
    for i in range(len(actions)):
        action_name = actions[i]
        is_equilibrium = True
        
        # Payoff if both players choose action i
        equilibrium_payoff = utility_matrix_alice[i, i]
        
        print(f"--- Checking strategy ({action_name}, {action_name}) ---")
        print(f"If both play {action_name}, Alice's payoff is {equilibrium_payoff}.")
        print("Now, let's see if Alice has an incentive to deviate, assuming Bob plays " + action_name + ".")
        
        # Check if Alice can get a better payoff by deviating
        for j in range(len(actions)):
            if i == j:
                continue
            
            deviating_payoff = utility_matrix_alice[j, i]
            print(f"  - If Alice deviates to {actions[j]}, her payoff would be {deviating_payoff}.")
            if deviating_payoff > equilibrium_payoff:
                print(f"    Since {deviating_payoff} > {equilibrium_payoff}, deviating is profitable. ({action_name}, {action_name}) is NOT a Nash Equilibrium.")
                is_equilibrium = False
                break
        
        if is_equilibrium:
            print(f"  No deviation is profitable for Alice. Therefore, ({action_name}, {action_name}) is a symmetric Nash Equilibrium.")
            symmetric_equilibria.append({'action': action_name, 'payoff': equilibrium_payoff})
        print("-" * 35 + "\n")

    print("Conclusion:")
    if not symmetric_equilibria:
        print("There are no symmetric pure-strategy Nash Equilibria.")
        # Note: A deeper analysis showed no symmetric mixed-strategy NE either, but that's beyond this code's scope.
    elif len(symmetric_equilibria) == 1:
        best_equilibrium = symmetric_equilibria[0]
        print(f"The only symmetric Nash Equilibrium is ({best_equilibrium['action']}, {best_equilibrium['action']}).")
        print(f"Superrational players will converge on this single best symmetric choice.")
        print(f"\nAlice's expected utility is her payoff in this equilibrium.")
        
        final_payoff = best_equilibrium['payoff']
        
        # As requested, printing the numbers in the final equation.
        # The equation for utility in the (Rest, Rest) state is direct.
        print("\nFinal Equation:")
        print(f"Alice's Utility = Payoff(Alice Rests, Bob Rests) = {final_payoff}")

    else:
        # If there were multiple symmetric NEs, they'd pick the one with the highest payoff.
        best_equilibrium = max(symmetric_equilibria, key=lambda x: x['payoff'])
        print(f"There are multiple symmetric Nash Equilibria. Superrational players choose the one with the highest payoff.")
        print(f"The best symmetric Nash Equilibrium is ({best_equilibrium['action']}, {best_equilibrium['action']}) with a payoff of {best_equilibrium['payoff']}.")
        
        final_payoff = best_equilibrium['payoff']
        print("\nFinal Equation:")
        print(f"Alice's Utility = {final_payoff}")


solve_utility()
<<<0>>>