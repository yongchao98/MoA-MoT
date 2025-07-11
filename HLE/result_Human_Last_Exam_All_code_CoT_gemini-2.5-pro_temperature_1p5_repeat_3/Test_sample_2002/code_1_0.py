def solve_token_game():
    """
    This function determines and prints the values of T for which the second
    player has a winning strategy in the described token game.

    The second player (P2) has a winning strategy if and only if the
    initial number of tokens, T, is a Fibonacci number. This is a known
    result from the analysis of a game called Fibonacci Nim, which matches
    the rules described. The strategy is based on Zeckendorf's theorem.

    This code generates the Fibonacci numbers up to a given limit, which
    are the values of T for which P2 wins.
    """
    limit = 1000
    print(f"The second player has a winning strategy if and only if T is a Fibonacci number.")
    print(f"The values of T for which the second player wins (up to {limit}) are:")

    # The Fibonacci sequence for Zeckendorf's representation, and thus for
    # this game's solution, is typically defined as F_1=1, F_2=2, F_3=3, ...
    a, b = 1, 2
    winning_t_values = []
    while a <= limit:
        winning_t_values.append(a)
        a, b = b, a + b

    # There isn't a single equation, so we print the list of numbers.
    # Each number in this list represents a value of T where P2 wins.
    print(", ".join(map(str, winning_t_values)))

solve_token_game()