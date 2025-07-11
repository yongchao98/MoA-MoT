import math

def check_for_psne(n_players):
    """
    Checks for Pure Strategy Nash Equilibria (PSNE) in a specific
    n-player symmetric game with 2 actions.

    The game is defined by the payoff differences f(k), which represents the
    incentive to switch to action '1' when k other players choose '1'.
    """
    print(f"Analyzing a {n_players}-player symmetric game...")
    print("-" * 30)

    # Let's define a game where the incentive to choose action '1' is always negative.
    # f(k) = u(1, k) - u(0, k)
    # Let's set f(k) = -1 for all k.
    # f is an array of length n_players, representing f(0), f(1), ..., f(n-1).
    f = [-1] * n_players

    total_psne_profiles = 0
    
    # Iterate through all possible numbers of players 'm' choosing action '1'.
    for m in range(n_players + 1):
        is_psne = False
        print(f"\nChecking profile type where m={m} players choose action '1':")

        # Case 1: All players choose action '0' (m=0)
        if m == 0:
            # Condition: a player switching to '1' would not be better off.
            # They would see k=0 others playing '1'.
            # The condition is u(0,0) >= u(1,0), which means f(0) <= 0.
            k = 0
            val = f[k]
            print(f"Checking condition for m=0: f[0] <= 0")
            print(f"We have f[0] = {val}. The condition {-1 <= 0} is {val <= 0}.")
            if val <= 0:
                is_psne = True
        
        # Case 2: All players choose action '1' (m=n)
        elif m == n_players:
            # Condition: a player switching to '0' would not be better off.
            # They would see k=n-1 others playing '1'.
            # The condition is u(1,n-1) >= u(0,n-1), which means f(n-1) >= 0.
            k = n_players - 1
            val = f[k]
            print(f"Checking condition for m={n_players}: f[{k}] >= 0")
            print(f"We have f[{k}] = {val}. The condition {val >= 0} is {val >= 0}.")
            if val >= 0:
                is_psne = True

        # Case 3: A mix of players (0 < m < n)
        else:
            # Condition 1: Players choosing '1' don't want to switch to '0'.
            # They see k=m-1 others playing '1'. Condition: f(m-1) >= 0.
            k1 = m - 1
            val1 = f[k1]
            print(f"Checking condition 1 (no '1'->'0' switch): f[{k1}] >= 0")
            print(f"We have f[{k1}] = {val1}. The condition {val1 >= 0} is {val1 >= 0}.")
            
            # Condition 2: Players choosing '0' don't want to switch to '1'.
            # They see k=m others playing '1'. Condition: f(m) <= 0.
            k2 = m
            val2 = f[k2]
            print(f"Checking condition 2 (no '0'->'1' switch): f[{k2}] <= 0")
            print(f"We have f[{k2}] = {val2}. The condition {val2 <= 0} is {val2 <= 0}.")
            
            if val1 >= 0 and val2 <= 0:
                is_psne = True

        # If conditions are met, count the number of profiles of this type
        if is_psne:
            if m == 0 or m == n_players:
                num_profiles = 1
            else:
                num_profiles = math.comb(n_players, m)
            
            print(f"--> SUCCESS: Profile type m={m} is a PSNE.")
            print(f"    This corresponds to {num_profiles} strategy profile(s).")
            total_psne_profiles += num_profiles
        else:
            print(f"--> FAILURE: Profile type m={m} is not a PSNE.")

    print("\n" + "-" * 30)
    print(f"Total number of Pure Strategy Nash Equilibria found in this game: {total_psne_profiles}")

# Let's run the check for a 5-player game
check_for_psne(n_players=5)