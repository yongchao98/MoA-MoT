def solve_token_game():
    """
    Determines and prints the initial number of tokens, T, for which the
    second player has a winning strategy in the described game.

    The analysis shows that this is a case of Fibonacci Nim, where the
    second player wins if and only if the starting number of tokens T
    is a Fibonacci number.
    """

    print("The second player has a winning strategy if and only if the initial number of tokens, T, is a Fibonacci number.")
    print("The Fibonacci sequence is typically defined as F(n) = F(n-1) + F(n-2).")
    print("Using the sequence starting with 1, 1, 2, 3, 5, ... we find the winning T values.")

    limit = 1000  # We will list the winning values of T up to this limit.
    print(f"\nThe values of T up to {limit} for which the second player wins are:")

    fib_numbers = []
    a, b = 1, 1
    # Generate Fibonacci numbers up to the limit
    while a <= limit:
        # Add to list only if not a duplicate (to handle the first two 1s)
        if a not in fib_numbers:
            fib_numbers.append(a)
        a, b = b, a + b

    # Print each of the identified Fibonacci numbers
    for t_value in fib_numbers:
        print(t_value)

solve_token_game()