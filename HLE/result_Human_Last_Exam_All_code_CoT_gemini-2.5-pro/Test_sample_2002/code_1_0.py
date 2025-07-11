def find_winning_tokens_for_p2():
    """
    This function identifies and prints the initial number of tokens, T, for which the second player (P2) has a winning strategy.

    The solution to this impartial game is based on its connection to the Fibonacci sequence. The second player has a
    winning strategy if and only if the initial number of tokens, T, is a Fibonacci number.

    This script generates and prints the unique Fibonacci numbers up to a limit (1000), which correspond to the
    winning values of T for the second player.
    """
    print("The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print("The winning values up to 1000 are:")

    # The sequence of winning T values is 1, 2, 3, 5, 8, ...
    # We can generate this by starting our Fibonacci sequence with a=1 and b=2.
    a, b = 1, 2
    
    # Print the first Fibonacci number, as our loop starts with the second.
    print(a)
    
    # Generate and print subsequent Fibonacci numbers up to the limit.
    limit = 1000
    while b <= limit:
        print(b)
        # Standard Fibonacci update step
        a, b = b, a + b

find_winning_tokens_for_p2()