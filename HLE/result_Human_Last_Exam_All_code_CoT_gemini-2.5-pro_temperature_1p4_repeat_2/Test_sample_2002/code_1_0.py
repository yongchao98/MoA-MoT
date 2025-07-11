def find_winning_tokens_for_p2():
    """
    This function determines and prints the initial number of tokens T for which
    the second player has a winning strategy in the described game.
    """

    print("The second player has a winning strategy if and only if the initial number of tokens, T, is a Fibonacci number.")
    print("This is a known result from combinatorial game theory for a game called Fibonacci Nim.")
    print("\nA position is a 'losing' position (for the current player) if the number of tokens is a Fibonacci number.")
    print("If T is a Fibonacci number, any move Player 1 makes will leave a non-Fibonacci number of tokens. Player 2 can then always make a move to leave a Fibonacci number for Player 1, eventually winning the game.")
    
    # For this game, the sequence F_n = F_{n-1} + F_{n-2} starts with F_1=1, F_2=2.
    print("\nThe Fibonacci sequence relevant to this game starts with 1 and 2.")
    print("\nHere are the values of T (up to 1000) for which the second player wins:")

    fib_numbers = []
    a, b = 1, 2
    limit = 1000

    # The first number in the sequence
    if a <= limit:
        fib_numbers.append(a)

    # Generate subsequent Fibonacci numbers up to the limit
    while b <= limit:
        fib_numbers.append(b)
        # The equation for generating the next number is b = a + (b-a).
        # We update a and b for the next iteration.
        a, b = b, a + b

    # Print all the found numbers
    print(", ".join(map(str, fib_numbers)))

find_winning_tokens_for_p2()