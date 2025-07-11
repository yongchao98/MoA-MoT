def find_winning_t_for_second_player():
    """
    This function identifies and prints the initial number of tokens, T,
    for which the second player has a winning strategy in the given game.

    The winning strategy is based on the properties of Fibonacci Nim.
    The second player wins if and only if T is a Fibonacci number.
    """
    print("The second player has a winning strategy if the initial number of tokens, T, is a Fibonacci number.")
    print("\nThe values of T for which the second player wins are the numbers in the Fibonacci sequence.")
    print("Here are the winning values for T up to 200:")

    # We list the unique positive Fibonacci numbers. The sequence of pile sizes is 1, 2, 3, 5, ...
    limit = 200
    a, b = 0, 1
    fib_sequence = []
    
    # Generate Fibonacci numbers up to the limit
    while b <= limit:
        # The game requires a positive number of tokens.
        if b > 0:
            # We are interested in the distinct possible values for T.
            if b not in fib_sequence:
                fib_sequence.append(b)
        
        # Standard Fibonacci generation
        a, b = b, a + b

    # Print the resulting list of numbers
    # The instruction "output each number in the final equation" is interpreted
    # as printing the numbers that form the solution set.
    print(*fib_sequence, sep=", ")
    print("...and so on, following the Fibonacci sequence.")

find_winning_t_for_second_player()