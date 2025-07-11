def main():
    """
    This program determines for which initial number of tokens, T, the second player
    has a winning strategy.

    The game described is a variant of the game of Nim, often called Fibonacci Nim.
    In this game, a position is a losing position for the current player if and only
    if the number of tokens is a Fibonacci number.

    Since the first player starts the game with T tokens, they are in a losing position
    if T is a Fibonacci number. Therefore, the second player has a winning strategy
    if and only if T belongs to the set of Fibonacci numbers {1, 2, 3, 5, 8, ...}.

    This script will now generate and print these numbers.
    """

    print("The second player has a winning strategy if T is a Fibonacci number.")
    print("The values of T for which the second player wins are:")

    # We set an arbitrary limit to display the first few winning T values.
    # The pattern continues indefinitely for all Fibonacci numbers.
    limit = 1000

    # The sequence of losing positions (and thus winning T for Player 2)
    # is 1, 2, 3, 5, 8, ...
    # We can generate this by starting our Fibonacci sequence with a=1 and b=2.

    a, b = 1, 2

    # Print the first term in the sequence
    if a <= limit:
        print(a)

    # Generate and print the rest of the terms up to the limit
    while b <= limit:
        print(b)
        # Standard Fibonacci calculation for the next term
        a, b = b, a + b

if __name__ == "__main__":
    main()