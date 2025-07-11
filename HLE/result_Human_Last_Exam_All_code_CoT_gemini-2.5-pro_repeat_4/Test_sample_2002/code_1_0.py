def generate_winning_t_values(limit):
    """
    This function generates the values of T for which the second player has a winning strategy.
    These values are the Fibonacci numbers.

    The game rules are:
    - A pile of T tokens.
    - P1 takes k1 tokens, where 0 < k1 < T.
    - Subsequent players can take k_new tokens, where 0 < k_new <= 2 * k_prev.
    - The player who takes the last token wins.

    The second player has a winning strategy if and only if T is a Fibonacci number.
    The relevant Fibonacci sequence for this game is 1, 2, 3, 5, 8, ...
    """
    print("The second player has a winning strategy if T is a Fibonacci number.")
    print(f"The values of T for which the second player wins (up to {limit}) are:")

    # We use the game-theoretic Fibonacci sequence which starts with 1, 2.
    # This corresponds to F_2, F_3, ... in the standard sequence F_0=0, F_1=1.
    # T=1 is a winning position for P2 (P1 has no valid moves).
    # T=2 is a winning position for P2 (P1 must take 1, P2 takes the last 1).
    a, b = 1, 2
    if a <= limit:
        print(a)
    while b <= limit:
        print(b)
        a, b = b, a + b

# Generate the winning numbers for T up to 1000.
generate_winning_t_values(1000)