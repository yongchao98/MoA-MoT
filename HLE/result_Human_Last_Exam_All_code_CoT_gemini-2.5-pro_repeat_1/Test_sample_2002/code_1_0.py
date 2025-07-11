def find_p2_winning_values(limit):
    """
    This function finds the values of T for which the second player has a
    winning strategy in the described token game. These values are the
    Fibonacci numbers.

    Args:
        limit (int): The upper bound for the values of T to check.
    """
    print(f"The values of T up to {limit} for which the second player has a winning strategy are the Fibonacci numbers.")
    print("These numbers are:")
    
    # The Fibonacci sequence is typically 1, 1, 2, 3, 5, ...
    # The set of winning positions for P2 is {1, 2, 3, 5, ...}.
    # We can generate this set.
    a, b = 0, 1
    fib_numbers = []
    while True:
        # This generates the standard Fibonacci sequence
        a, b = b, a + b
        if a > limit:
            break
        # We add to a set to handle the initial duplicate 1
        if a not in fib_numbers:
            fib_numbers.append(a)

    for num in fib_numbers:
        print(num)

# You can change this limit to find more values of T
upper_limit = 200
find_p2_winning_values(upper_limit)