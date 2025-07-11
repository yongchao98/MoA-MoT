def find_spne():
    """
    Analyzes the provided Nash Equilibria to determine which are subgame perfect.
    """
    
    # Payoffs for Player 1 in the subgame starting after the history (D, u).
    # Player 1 chooses between F and B.
    payoffs_subgame = {
        'F': -1,
        'B': 0
    }
    
    print("Step 1: Identify the proper subgame.")
    print("A proper subgame starts at Player 1's decision node after the sequence of moves (D, u).")
    print("-" * 30)
    
    print("Step 2: Analyze the rational choice in the subgame.")
    optimal_choice_p1 = max(payoffs_subgame, key=payoffs_subgame.get)
    print(f"In this subgame, Player 1 chooses between F and B.")
    print(f"Payoff for F: {payoffs_subgame['F']}")
    print(f"Payoff for B: {payoffs_subgame['B']}")
    print(f"Since {payoffs_subgame['B']} > {payoffs_subgame['F']}, the rational choice for Player 1 is '{optimal_choice_p1}'.")
    print("-" * 30)

    # The list of pure strategy Nash Equilibria given in the problem.
    # Format: (Player1_Strategy, Player2_Strategy)
    # P1's strategy is a string: {Initial Move}{Move after U}{Move after (D,u)}
    nash_equilibria = {
        "A": ("U RF", "u"),
        "B": ("DLF", "d"), # Not a real option, just for mapping
        "C": ("DRF", "d"), # Not a real option, just for mapping
        "D": ("U RB", "u")
    }
    
    # Let's use the actual NEs from the problem statement
    nes_from_problem = [("U RF", "u"), ("U RB", "u"), ("DLF", "d"), ("DRF", "d")]
    
    spne_list = []
    
    print("Step 3: Check each Nash Equilibrium against the subgame perfection criterion.")
    print("An SPNE must specify the optimal action in every subgame.")
    print(f"Therefore, Player 1's strategy must specify '{optimal_choice_p1}' for the subgame.\n")
    
    for ne in nes_from_problem:
        p1_strategy = ne[0]
        # The last character of P1's strategy string is the move in the subgame.
        action_in_subgame = p1_strategy[-1]
        
        print(f"Checking NE: {ne}")
        print(f"  Player 1's strategy is '{p1_strategy}'.")
        print(f"  In the subgame, this strategy specifies the action '{action_in_subgame}'.")
        
        if action_in_subgame == optimal_choice_p1:
            print(f"  This action is optimal. Therefore, {ne} is a Subgame Perfect Nash Equilibrium.")
            spne_list.append(ne)
        else:
            print(f"  This action is not optimal (Player 1 would prefer '{optimal_choice_p1}' for a payoff of {payoffs_subgame[optimal_choice_p1]} instead of {payoffs_subgame[action_in_subgame]}).")
            print(f"  Therefore, {ne} is NOT a Subgame Perfect Nash Equilibrium.")
        print("-" * 20)
        
    print("\nFinal Conclusion:")
    print(f"The Subgame Perfect Nash Equilibrium from the given list is: {spne_list[0]}")

find_spne()