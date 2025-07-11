def solve_token_game():
    """
    This function finds the values of T for which the second player has a winning strategy.

    The winning strategy for the second player exists if and only if the initial number
    of tokens, T, is a Fibonacci number. This code generates Fibonacci numbers
    up to a specified limit.
    """
    limit = 1000
    
    print(f"The second player has a winning strategy if and only if T is a Fibonacci number.")
    print(f"Here are the values of T up to {limit}:")

    # Generate Fibonacci numbers. The game is usually considered for T>=2.
    # The sequence of losing positions starts with F_3=2, F_4=3, ... if we use F_1=1, F_2=1
    # We will generate the sequence 2, 3, 5, 8, ...
    a, b = 2, 3
    
    fib_numbers = []
    
    # Handle the first two Fibonacci numbers relevant to the game (T>=2)
    if 2 <= limit:
        fib_numbers.append(2)
    
    # Generate the rest of the sequence
    while b <= limit:
        fib_numbers.append(b)
        a, b = b, a + b
    
    # Print the result
    print(', '.join(map(str, fib_numbers)))

solve_token_game()