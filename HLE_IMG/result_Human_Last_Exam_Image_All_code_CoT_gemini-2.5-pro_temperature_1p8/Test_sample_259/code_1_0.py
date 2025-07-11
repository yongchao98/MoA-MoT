def find_spne():
    """
    This function identifies the Subgame Perfect Nash Equilibria (SPNE) from a given
    list of Nash Equilibria (NE) by applying the principle of backward induction.
    """

    # Payoffs are represented as (Player 1's payoff, Player 2's payoff)

    # Step 1: Analyze the only proper subgame.
    # The subgame starts after Player 1 chooses 'D' and Player 2 chooses 'u'.
    # In this subgame, Player 1 chooses between 'F' and 'B'.
    payoff_F = (-1, -1)
    payoff_B = (0, 4)
    p1_payoff_F = payoff_F[0]
    p1_payoff_B = payoff_B[0]

    print("Analyzing the subgame where Player 1 chooses between F and B:")
    print(f"If Player 1 chooses F, their payoff is {p1_payoff_F}.")
    print(f"If Player 1 chooses B, their payoff is {p1_payoff_B}.")

    # Determine Player 1's optimal choice in this subgame.
    if p1_payoff_B > p1_payoff_F:
        optimal_choice_p1 = 'B'
        print(f"Since {p1_payoff_B} > {p1_payoff_F}, the rational choice for Player 1 in this subgame is '{optimal_choice_p1}'.")
    else:
        optimal_choice_p1 = 'F'
        print(f"Since {p1_payoff_F} >= {p1_payoff_B}, the rational choice for Player 1 in this subgame is '{optimal_choice_p1}'.")

    # Step 2: Filter the given Nash Equilibria.
    # An NE is an SPNE only if it specifies the optimal action in every subgame.
    # Player 1's strategy is a tuple: (Action at root, Action after U, Action after D->u)
    # The third element must be the optimal choice 'B'.
    
    nash_equilibria = {
        "(U RF, u)": ('U', 'R', 'F'),
        "(U RB, u)": ('U', 'R', 'B'),
        "(DLF, d)": ('D', 'L', 'F'),
        "(DRF, d)": ('D', 'R', 'F')
    }

    print("\nChecking the provided Nash Equilibria:")
    spne = []
    for ne_str, p1_strategy in nash_equilibria.items():
        # The choice in the subgame is the third element of Player 1's strategy tuple.
        choice_in_subgame = p1_strategy[2]
        if choice_in_subgame == optimal_choice_p1:
            print(f"- {ne_str}: Player 1's strategy is {p1_strategy}. The choice in the subgame ('{choice_in_subgame}') is optimal. This IS a Subgame Perfect Nash Equilibrium.")
            spne.append(ne_str)
        else:
            print(f"- {ne_str}: Player 1's strategy is {p1_strategy}. The choice in the subgame ('{choice_in_subgame}') is not optimal. This is NOT a Subgame Perfect Nash Equilibrium.")

    print("\nThe Subgame Perfect Nash Equilibrium (or equilibria) in pure strategies is:")
    for s in spne:
        print(s)

find_spne()