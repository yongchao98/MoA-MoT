def find_winning_T_for_player2():
    """
    This script calculates and prints the initial number of tokens, T, for which 
    the second player has a winning strategy in the game of Fibonacci Nim.

    The second player wins if and only if the starting number of tokens T
    is a Fibonacci number. The game requires T > 1 to start.
    """

    print("The second player has a winning strategy if T is a Fibonacci number greater than 1.")
    print("This is an infinite set of numbers. The first few values of T are:")

    limit = 200  # We will generate numbers up to this limit as examples.
    
    # We start with the Fibonacci sequence definition F(n) = F(n-1) + F(n-2)
    # Using the standard seeds F(1)=1, F(2)=1, the next numbers are 2, 3, 5, ...
    a, b = 1, 1
    
    winning_T = []
    
    while True:
        # Calculate the next Fibonacci number
        next_fib = a + b
        
        # Stop if we exceed the limit
        if next_fib > limit:
            break
            
        winning_T.append(next_fib)
        
        # Update the sequence
        a, b = b, next_fib
        
    # Print the resulting list of numbers
    print(', '.join(map(str, winning_T)))

find_winning_T_for_player2()