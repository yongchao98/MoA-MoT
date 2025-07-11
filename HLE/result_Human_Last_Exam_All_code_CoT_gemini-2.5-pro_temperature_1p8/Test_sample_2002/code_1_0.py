def find_p2_winning_T(limit):
    """
    This function finds and prints the values of T for which the second player has a
    winning strategy. These values are Fibonacci numbers.

    The problem states rules for an impartial game. The positions from which
    the second player has a winning strategy are the "cold" positions. For this
    specific game, these correspond to the Fibonacci numbers.

    We use the standard Fibonacci sequence definition: F_0=0, F_1=1, F_2=1, ...
    The winning T values for P2 are F_n for n >= 3 (i.e., T >= 2).
    """
    print(f"The values of T up to {limit} for which the second player has a winning strategy are:")
    # Initialize with the first two Fibonacci numbers
    a, b = 0, 1
    while True:
        # Calculate the next Fibonacci number
        fib = a + b
        if fib > limit:
            break
        # P2 wins for T that are Fibonacci numbers (starting from T=2)
        if fib >= 2:
            print(fib)
        a = b
        b = fib

# Set a limit for T
upper_limit = 200
find_p2_winning_T(upper_limit)