def find_winning_T_for_player2():
    """
    This function identifies the initial number of tokens T for which the
    second player has a winning strategy in the described token game.

    The problem is a variant of the game of Nim known as Fibonacci Nim.
    The second player has a winning strategy if and only if the initial
    number of tokens, T, is a Fibonacci number.

    The sequence of Fibonacci numbers relevant to this game is:
    1, 2, 3, 5, 8, 13, 21, ...
    where each number is the sum of the two preceding ones.

    This script will generate and print these numbers up to a practical limit (e.g., 1000).
    """
    print("The second player has a winning strategy if and only if the initial number of tokens T is a Fibonacci number.")
    print("The values of T are members of the Fibonacci sequence (1, 2, 3, 5, 8, ...).")
    print("\nHere are the winning values for T up to 1000:")

    # We use the Fibonacci sequence starting with 1, 2.
    a, b = 1, 2
    fib_numbers = []
    limit = 1000

    while a <= limit:
        fib_numbers.append(a)
        a, b = b, a + b
    
    # Printing the numbers as requested.
    # The prompt asked to "output each number in the final equation", which is
    # interpreted here as printing the elements of the resulting set of numbers.
    print(*fib_numbers, sep=", ")

find_winning_T_for_player2()