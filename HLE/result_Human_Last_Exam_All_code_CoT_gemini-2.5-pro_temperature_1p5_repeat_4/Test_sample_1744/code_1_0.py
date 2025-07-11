import math

def count_ne_profiles(n, f, game_name):
    """
    Calculates the number of pure strategy Nash Equilibrium profiles for a given game.
    
    Args:
        n (int): The number of players.
        f (list): A list of n floats representing the function f(k) = u(1,k) - u(0,k)
                  for k from 0 to n-1.
        game_name (str): The name of the game scenario for printing.
    """
    print(f"--- Analyzing Game: {game_name} for n={n} ---")
    total_ne_profiles = 0
    
    # Check for NE where K=0 players play action 1 (all play 0)
    # Condition: f(0) <= 0
    k_val = 0
    print(f"Checking for NE type K={k_val} (all players play action 0):")
    print(f"  Condition: f[0] <= 0")
    print(f"  Check: {f[0]} <= 0 is {f[0] <= 0}")
    if f[0] <= 0:
        num_profiles = math.comb(n, k_val)
        print(f"  Result: K=0 is a NE type. Found {num_profiles} profile.")
        total_ne_profiles += num_profiles
    else:
        print("  Result: K=0 is not a NE type.")

    # Check for NE where 0 < K < n players play action 1
    for K in range(1, n):
        print(f"\nChecking for NE type K={K} ({K} players play 1, {n-K} players play 0):")
        # Condition: f(K-1) >= 0 AND f(K) <= 0
        print(f"  Condition: f[{K-1}] >= 0 AND f[{K}] <= 0")
        print(f"  Check: {f[K-1]} >= 0 is {f[K-1] >= 0}, and {f[K]} <= 0 is {f[K] <= 0}")
        if f[K-1] >= 0 and f[K] <= 0:
            num_profiles = math.comb(n, K)
            print(f"  Result: K={K} is a NE type. Found C({n}, {K}) = {num_profiles} profiles.")
            total_ne_profiles += num_profiles
        else:
            print(f"  Result: K={K} is not a NE type.")
            
    # Check for NE where K=n players play action 1 (all play 1)
    # Condition: f(n-1) >= 0
    k_val = n
    print(f"\nChecking for NE type K={k_val} (all players play action 1):")
    print(f"  Condition: f[{n-1}] >= 0")
    print(f"  Check: {f[n-1]} >= 0 is {f[n-1] >= 0}")
    if f[n-1] >= 0:
        num_profiles = math.comb(n, k_val)
        print(f"  Result: K={n} is a NE type. Found {num_profiles} profile.")
        total_ne_profiles += num_profiles
    else:
        print("  Result: K={n} is not a NE type.")
    
    print(f"\n--- Total number of PSNE profiles for '{game_name}': {total_ne_profiles} ---\n")
    return total_ne_profiles

def main():
    n_players = 5

    # Case 1: Game where action 0 is always preferred (e.g., a coordination game on action 0)
    # We expect exactly 1 NE profile: (0, 0, 0, 0, 0)
    f_case1 = [-1.0] * n_players  # f(k) = -1 for all k
    count_ne_profiles(n_players, f_case1, "Action 0 is Dominant")

    # Case 2: Game where action 1 is always preferred
    # We expect exactly 1 NE profile: (1, 1, 1, 1, 1)
    f_case2 = [1.0] * n_players  # f(k) = 1 for all k
    count_ne_profiles(n_players, f_case2, "Action 1 is Dominant")

    # Case 3: Anti-coordination game (e.g., Hawk-Dove style)
    # Let f(k) increase with k, e.g., f(k) = k - (n-1)/2
    # f = [-2, -1, 0, 1, 2] for n=5
    # We expect NE at K=2 (f[1]=-1<0, f[2]=0<=0) -- no
    # Let's adjust, NE condition: f[K-1]>=0 and f[K]<=0
    # Let's have f be decreasing. f = [2, 1, 0.5, -0.5, -1]
    # We should have a NE at K=3. f[2]=0.5>=0 and f[3]=-0.5<=0
    f_case3 = [2.0, 1.0, 0.5, -0.5, -1.0]
    count_ne_profiles(n_players, f_case3, "Anti-Coordination")
    
if __name__ == "__main__":
    main()