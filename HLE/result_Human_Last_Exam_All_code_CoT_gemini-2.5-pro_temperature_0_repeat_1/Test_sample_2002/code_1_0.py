def find_winning_t():
    """
    This function identifies and prints the initial number of tokens, T,
    for which the second player has a winning strategy in the described game.

    The game is a well-known impartial game called Fibonacci Nim. The second
    player has a winning strategy if and only if T is a Fibonacci number.
    The relevant Fibonacci sequence starts with 1 and 2.
    """
    limit = 1000  # We will find the winning values for T up to this limit.

    print(f"The second player has a winning strategy if T is a Fibonacci number.")
    print(f"The first few winning values for T up to {limit} are:")

    a, b = 1, 2

    # Print the first two Fibonacci numbers in the sequence
    if a <= limit:
        print(a)
    if b <= limit:
        print(b)

    # Generate and print subsequent Fibonacci numbers using the recurrence relation
    while True:
        next_fib = a + b
        if next_fib > limit:
            break
        print(f"{next_fib} = {b} + {a}")
        a, b = b, next_fib

find_winning_t()