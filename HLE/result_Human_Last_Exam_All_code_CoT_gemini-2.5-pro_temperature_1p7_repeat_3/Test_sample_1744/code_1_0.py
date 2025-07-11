import math

def nCr(n, r):
    """Calculate the number of combinations (n choose r)."""
    # This function is used to count how many strategy profiles correspond to a given 'k'.
    # For example, if n=5 and k=2, there are C(5,2)=10 profiles with two 'A's and three 'B's.
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

def find_psne_profiles(n, d_func):
    """
    Finds the number of pure strategy Nash equilibria (PSNE) profiles for a given
    n-player symmetric game with 2 actions.
    
    Args:
        n (int): The number of players.
        d_func (function): A function d(m) = u(A, m) - u(B, m), representing the
                           incentive to choose 'A' when m others choose 'A'.
    """
    total_psne_profiles = 0
    print(f"Analyzing a symmetric game with n={n} players and actions 'A' and 'B'.")
    print("A pure strategy profile is a Nash Equilibrium if no player can improve their payoff by unilaterally changing their action.")
    print("We check for each possible number of 'A' players, k from 0 to n.\n")
    
    # Case k=0: All players choose 'B'
    k = 0
    # The condition for this to be a PSNE is that a single player switching to 'A'
    # would not be better off. This player sees m=0 others playing 'A'.
    # Condition: u(B, 0) >= u(A, 0)  --->  d(0) <= 0
    d_0 = d_func(0)
    is_psne_k0 = (d_0 <= 0)
    print(f"Checking for k={k} (all 'B'):")
    print(f"  Condition: d(0) <= 0. We have d(0) = {d_0}. The condition is {is_psne_k0}.")
    if is_psne_k0:
        num_profiles = nCr(n, k)
        total_psne_profiles += num_profiles
        print(f"  Result: This is a PSNE type. Found {num_profiles} profile.")
    else:
        print("  Result: Not a PSNE type.")
    print("-" * 20)

    # Case 0 < k < n: Mixed population of 'A' and 'B' players
    for k in range(1, n):
        # Condition 1: An 'A' player does not want to switch to 'B'. They see k-1 other 'A' players.
        # Condition: u(A, k-1) >= u(B, k-1)  ---> d(k-1) >= 0
        d_k_minus_1 = d_func(k - 1)
        cond1 = (d_k_minus_1 >= 0)
        
        # Condition 2: A 'B' player does not want to switch to 'A'. They see k other 'A' players.
        # Condition: u(B, k) >= u(A, k) ---> d(k) <= 0
        d_k = d_func(k)
        cond2 = (d_k <= 0)
        
        print(f"Checking for k={k} ({k} 'A', {n-k} 'B'):")
        print(f"  Condition (A-player): d({k-1}) >= 0. We have d({k-1}) = {d_k_minus_1}. Condition is {cond1}.")
        print(f"  Condition (B-player): d({k}) <= 0. We have d({k}) = {d_k}. Condition is {cond2}.")
        
        is_psne_k = cond1 and cond2
        if is_psne_k:
            num_profiles = nCr(n, k)
            total_psne_profiles += num_profiles
            print(f"  Result: Both conditions met. This is a PSNE type. Found {num_profiles} profiles.")
        else:
            print("  Result: Not a PSNE type.")
        print("-" * 20)

    # Case k=n: All players choose 'A'
    k = n
    # The condition is that a player switching to 'B' is not better off.
    # This player sees m=n-1 others playing 'A'.
    # Condition: u(A, n-1) >= u(B, n-1) ---> d(n-1) >= 0
    d_n_minus_1 = d_func(n - 1)
    is_psne_kn = (d_n_minus_1 >= 0)
    print(f"Checking for k={n} (all 'A'):")
    print(f"  Condition: d({n-1}) >= 0. We have d({n-1}) = {d_n_minus_1}. The condition is {is_psne_kn}.")
    if is_psne_kn:
        num_profiles = nCr(n, k)
        total_psne_profiles += num_profiles
        print(f"  Result: This is a PSNE type. Found {num_profiles} profile.")
    else:
        print("  Result: Not a PSNE type.")
    print("-" * 20)
    
    print("\n" + "="*40)
    print(f"Total number of PSNE profiles found: {total_psne_profiles}")
    print("="*40)

def main():
    # Number of players
    n = 5
    
    # Define a game where 'B' is a strictly dominant strategy.
    # For this game, u(B, m) > u(A, m) for all m.
    # This implies d(m) = u(A, m) - u(B, m) is always negative.
    def dominant_strategy_game(m):
        return -1
        
    print("Demonstration with a specific game where the minimum is achieved.")
    find_psne_profiles(n, dominant_strategy_game)

if __name__ == '__main__':
    main()