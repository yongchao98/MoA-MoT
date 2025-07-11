import re

def find_spne():
    """
    This function analyzes the given game tree to find the Subgame Perfect Nash Equilibria (SPNE)
    from a provided list of Nash Equilibria (NE).
    """

    print("Step 1: Identifying Proper Subgames")
    print("A proper subgame starts at a decision node that is a singleton information set.")
    print("In the given game tree, the only proper subgame starts after the sequence of moves (D, u).")
    print("At this point, Player 1 must make a decision.\n")

    print("Step 2: Solving the Subgame using Backward Induction")
    print("In the subgame, Player 1 chooses between 'F' and 'B'.")
    
    # Payoffs for Player 1 in the subgame
    p1_payoff_F = -1
    p1_payoff_B = 0
    
    print(f"Player 1's payoff for choosing 'F' is: {p1_payoff_F}")
    print(f"Player 1's payoff for choosing 'B' is: {p1_payoff_B}")
    
    # Determine the rational choice for Player 1
    if p1_payoff_B > p1_payoff_F:
        rational_choice = 'B'
        print(f"Since {p1_payoff_B} > {p1_payoff_F}, a rational Player 1 must choose 'B' in this subgame.")
    else:
        rational_choice = 'F'
        print(f"Since {p1_payoff_F} >= {p1_payoff_B}, a rational Player 1 must choose 'F' in this subgame.")
    print("Therefore, for a strategy to be a Subgame Perfect Nash Equilibrium, it must specify that Player 1 chooses 'B'.\n")
    
    print("Step 3: Checking the Provided List of Nash Equilibria")
    # The given list of pure strategy Nash Equilibria
    # Player 1's strategy is represented as InitialMove-UpperMove-LowerMove, e.g., 'U RF'
    nes_pure = ["(U RF, u)", "(U RB, u)", "(DLF, d)", "(DRF, d)"]
    print(f"The list of NEs to check is: {nes_pure}\n")
    
    spne_list = []
    
    for ne in nes_pure:
        # Extract Player 1's strategy, which is the part before the comma, like "U RF"
        p1_strategy_full = ne.split(',')[0].strip()[1:] # e.g., 'U RF'
        
        # The action in the subgame is the last character of the second part of P1's strategy
        p1_action_in_subgame = p1_strategy_full.split(' ')[1][1]
        
        print(f"Analyzing NE: {ne}")
        print(f"Player 1's action in the subgame is '{p1_action_in_subgame}'.")
        
        if p1_action_in_subgame == rational_choice:
            print(f"This matches the rational choice '{rational_choice}'. This NE is Subgame Perfect.")
            spne_list.append(ne)
        else:
            print(f"This does not match the rational choice '{rational_choice}'. This NE is NOT Subgame Perfect.")
        print("-" * 20)

    print("\nStep 4: Final Conclusion")
    if not spne_list:
        print("No Subgame Perfect Nash Equilibria were found in the list.")
    else:
        print("The Subgame Perfect Nash Equilibrium (or Equilibria) is/are:")
        for spne in spne_list:
            print(spne)

find_spne()