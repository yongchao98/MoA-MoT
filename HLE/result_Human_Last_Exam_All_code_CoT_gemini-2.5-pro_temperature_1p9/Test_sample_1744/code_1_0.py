def find_minimum_psne_example(n=5):
    """
    This function analyzes a specific n-player, 2-action symmetric game
    to find the number of pure strategy Nash equilibria (PSNE).

    The game is an n-player "Prisoner's Dilemma" where action 'A' is always
    strictly dominant over action 'B'. This means a player always has an
    incentive to choose 'A', regardless of others' actions.

    We define an incentive function f(k) = Payoff(A, k) - Payoff(B, k),
    where k is the number of *other* players choosing A.
    For this game, we set f(k) = 1 for all k, indicating A is always better.
    """
    print(f"Analyzing a {n}-player game where action A is always preferred to B.")
    print("-" * 60)

    # f_values[k] corresponds to f(k)
    # We define a simple game where f(k) is always 1.
    f_values = [1] * n
    
    psne_list = []
    
    # Check for PSNE where k players choose action A.
    for k in range(n + 1):
        is_psne = False
        print(f"Checking if k={k} ( {k} players choose A, {n-k} choose B) is a PSNE:")
        
        # Case 1: k=0 (all players choose B)
        # Condition: f(0) <= 0
        if k == 0:
            f_0 = f_values[0]
            condition = (f_0 <= 0)
            print(f"  Condition: Incentive to switch to A, f(0), must be <= 0.")
            print(f"  Calculation: f(0) = {f_0}. Is {f_0} <= 0? {condition}")
            if condition:
                is_psne = True

        # Case 2: k=n (all players choose A)
        # Condition: f(n-1) >= 0
        elif k == n:
            f_n_minus_1 = f_values[n-1]
            condition = (f_n_minus_1 >= 0)
            print(f"  Condition: Incentive to stay with A, f(n-1), must be >= 0.")
            print(f"  Calculation: f({n-1}) = {f_n_minus_1}. Is {f_n_minus_1} >= 0? {condition}")
            if condition:
                is_psne = True

        # Case 3: 0 < k < n
        # Conditions: f(k-1) >= 0 AND f(k) <= 0
        else:
            f_k_minus_1 = f_values[k-1]
            f_k = f_values[k]
            condition1 = (f_k_minus_1 >= 0)
            condition2 = (f_k <= 0)
            print(f"  Condition 1 (A-players don't switch): f({k-1}) >= 0.")
            print(f"    Calculation: f({k-1}) = {f_k_minus_1}. Is {f_k_minus_1} >= 0? {condition1}")
            print(f"  Condition 2 (B-players don't switch): f({k}) <= 0.")
            print(f"    Calculation: f({k}) = {f_k}. Is {f_k} <= 0? {condition2}")
            if condition1 and condition2:
                is_psne = True

        if is_psne:
            psne_list.append(k)
            print("  Result: This is a PSNE.\n")
        else:
            print("  Result: This is not a PSNE.\n")

    print("-" * 60)
    print(f"Found {len(psne_list)} PSNE(s) for this game.")
    if len(psne_list) > 0:
      print(f"The equilibrium profile(s) (by number of players choosing A) are: {psne_list}")

find_minimum_psne_example()