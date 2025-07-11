def find_spne():
    """
    This function identifies the Subgame Perfect Nash Equilibria (SPNE)
    from a given list of Nash Equilibria (NE) for the specified game.
    """
    # The list of pure strategy Nash Equilibria given in the problem.
    # The format is (Player1_Strategy, Player2_Strategy).
    # Player 1's strategy is a string of 3 characters:
    # 1. Initial move (U/D)
    # 2. Move in the upper information set (L/R)
    # 3. Move in the lower subgame (F/B)
    nash_equilibria = [
        ("U RF", "u"),
        ("U RB", "u"),
        ("DLF", "d"),
        ("DRF", "d")
    ]

    print("Analyzing the subgame starting after moves (D, u):")
    p1_payoff_F = -1
    p1_payoff_B = 0
    print(f"Player 1's payoff for choosing F: {p1_payoff_F}")
    print(f"Player 1's payoff for choosing B: {p1_payoff_B}")
    
    if p1_payoff_B > p1_payoff_F:
        optimal_choice_in_subgame = 'B'
        print(f"Since {p1_payoff_B} > {p1_payoff_F}, the rational choice for Player 1 in this subgame is 'B'.")
    else:
        optimal_choice_in_subgame = 'F'
        print(f"Since {p1_payoff_F} >= {p1_payoff_B}, the rational choice for Player 1 in this subgame is 'F'.")

    print("\nAn SPNE must have Player 1 making the optimal choice in every subgame.")
    print(f"Filtering the list of Nash Equilibria for strategies where Player 1 chooses '{optimal_choice_in_subgame}' in the subgame.")
    
    subgame_perfect_equilibria = []
    for p1_strategy, p2_strategy in nash_equilibria:
        # The choice in the subgame is the third character of Player 1's strategy string.
        choice_in_subgame = p1_strategy[2]
        if choice_in_subgame == optimal_choice_in_subgame:
            subgame_perfect_equilibria.append(f"({p1_strategy}, {p2_strategy})")

    print("\nThe Subgame Perfect Nash Equilibrium (or Equilibria) in pure strategies is:")
    for spne in subgame_perfect_equilibria:
        print(spne)

find_spne()