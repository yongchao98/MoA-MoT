def find_winning_tokens_for_p2():
    """
    This function identifies and prints the number of tokens T for which
    the second player has a winning strategy in the described game.
    
    The winning strategy for the second player exists if and only if T
    is a Fibonacci number. This program generates Fibonacci numbers up to a
    predefined limit.
    """
    
    limit = 1000
    
    print("The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print(f"The winning values for T up to {limit} are:")
    
    # We use a set to store the Fibonacci numbers to automatically handle
    # the duplicate '1' at the beginning of the sequence (F_1=1, F_2=1).
    fib_numbers = set()
    
    # Initialize with the first two Fibonacci numbers
    a, b = 1, 1
    
    # Generate Fibonacci numbers up to the limit
    while a <= limit:
        fib_numbers.add(a)
        a, b = b, a + b
        
    # Print the unique, sorted Fibonacci numbers
    # The * operator unpacks the list for the print function.
    print(*sorted(list(fib_numbers)), sep=", ")

find_winning_tokens_for_p2()