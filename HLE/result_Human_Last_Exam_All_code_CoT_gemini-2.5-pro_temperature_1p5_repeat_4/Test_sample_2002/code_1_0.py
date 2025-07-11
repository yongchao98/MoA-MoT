def find_winning_tokens_for_p2(limit):
    """
    This function identifies the initial number of tokens T for which the second
    player has a winning strategy in the described token game.

    The game is a classic impartial game called Fibonacci Nim. The winning strategy
    is determined by the Fibonacci sequence. The second player has a winning
    strategy if and only if the initial number of tokens, T, is a Fibonacci number.

    This script generates Fibonacci numbers up to a specified limit to show
    which values of T are winning for the second player.
    """
    # We will generate Fibonacci numbers up to the given limit.
    # The sequence starts with 1, 2, 3, 5... The standard generation starting
    # with a=0, b=1 will produce the set {1, 2, 3, 5, ...}.
    fib_numbers = []
    a, b = 0, 1
    while b <= limit:
        # We add b to our list of winning numbers for player 2.
        # We use a set temporarily to store unique numbers and then convert to a sorted list
        # because the first two numbers in the sequence 1, 1, 2, ... are the same.
        if b > 0:
            fib_numbers.append(b)
        a, b = b, a + b
    
    # Remove duplicates and sort the list
    unique_fib_numbers = sorted(list(set(fib_numbers)))
    
    print("The second player has a winning strategy if and only if T is a Fibonacci number.")
    print(f"For a total number of tokens T up to {limit}, these values are:")
    
    # The instruction "output each number in the final equation" is interpreted
    # as printing the list of winning numbers for player 2.
    print(*unique_fib_numbers)

# Let's find the winning T values for the second player up to a limit of 200.
find_winning_tokens_for_p2(200)
