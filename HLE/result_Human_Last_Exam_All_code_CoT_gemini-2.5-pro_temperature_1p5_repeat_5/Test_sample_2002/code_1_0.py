def find_winning_tokens_for_p2(limit):
    """
    This function finds and prints the number of tokens T for which the second
    player has a winning strategy. These numbers are the Fibonacci numbers.

    Args:
    limit (int): The upper bound for the number of tokens T.
    """
    print(f"The values of T for which the second player has a winning strategy are the Fibonacci numbers.")
    print(f"Here are the Fibonacci numbers up to {limit}:")
    
    # We use the Fibonacci sequence starting with 1, 2, 3, 5...
    # as T=1 is also a winning position for P2.
    # We'll generate using a=1, b=1 to get the standard sequence 1, 1, 2, 3, ...
    # and handle the duplicate 1.
    a, b = 1, 1
    fib_numbers = []
    while a <= limit:
        if a not in fib_numbers:
             fib_numbers.append(a)
        # Using a temporary variable to calculate the next fibonacci number
        # is a more intuitive way to show the progression.
        next_fib = a + b
        a = b
        b = next_fib
        
    # We directly print the numbers as requested.
    output = ""
    for i, num in enumerate(fib_numbers):
        if i > 0:
            output += ", "
        output += str(num)
    print(output)

# Let's find the numbers for T up to 1000.
find_winning_tokens_for_p2(1000)