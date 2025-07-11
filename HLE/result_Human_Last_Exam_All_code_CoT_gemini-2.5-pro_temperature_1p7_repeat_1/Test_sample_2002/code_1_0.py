def find_winning_t_for_p2(limit):
    """
    This function finds the values of T for which the second player has a winning strategy.
    These values are the Fibonacci numbers.
    """
    
    # Using the standard Fibonacci sequence definition F_0=0, F_1=1, F_2=1, ...
    # The winning numbers for player 2 are 1, 2, 3, 5, 8, ...
    # which are all the Fibonacci numbers > 0.
    a, b = 0, 1
    fib_numbers = []
    
    while True:
        # Generate next Fibonacci number
        a, b = b, a + b
        if a > limit:
            break
        # The first number generated is 1. We don't want to add it twice.
        if a not in fib_numbers:
             fib_numbers.append(a)

    print(f"The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print(f"Here are the values up to {limit}:")
    # Using map to convert all integers in the list to string
    print(', '.join(map(str, fib_numbers)))

# We can set a reasonable limit for T.
find_winning_t_for_p2(500)
